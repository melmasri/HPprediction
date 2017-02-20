#!/usr/bin/env Rscript
library("optparse")

option_list = list(
    make_option(c("-f", "--file"), type="character",
                default=NULL, 
                help="dataset file name, or an identifiable portion of it, if multiple files the first one is used. The file should have the RData extension.",
                metavar="character"),
	make_option(c("-r", "--runtype"), type="character",
                default="10fold", 
                help="possible run types, ALL, 10fold, affinity, distonly, weighted [default= %default]",
                metavar="character"),
    make_option(c("-s", "--subset"), type="character",
                default=NULL, 
                help="possible subset choices, for GMP it is Carnivora, for and EID it is Rodentia.",
                metavar="character"),
    make_option("--output", type="character",
                default=NULL, 
                help="output directory, the default is constructed from time stamp and other inputs.",
                metavar="character"),
    make_option("--no.cycles", type="integer",
                default=20, 
                help="no. of cycles to run over the matrix [default= %default].",
                metavar="integer"),
    make_option("--no.cores", type="integer",
                default=5, 
                help="number of cores to use (1, 2, ..) [default= %default].",
                metavar="integer"),
    make_option("--icm.horiz", type="logical",
                default=TRUE, 
                help="block updates are done horizontally (ICM.HORIZ=TRUE), else diagonally [default= %default].",
                metavar="logical"),
    make_option("--email", type="character",
                default=NULL,
                help="end-of-run email notification is sent to the supplied email if mailx is installed and configured",
                metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

NO.CORES = opt$no.cores
ICM.HORIZ = opt$icm.horiz
SLICE = opt$no.cycles
## Setting up the file name 
if(any(grepl(opt$file, list.files(recursive=FALSE), ignore.case = TRUE))){
    DATAFILENAME  = grep(opt$file, list.files(recursive=FALSE, pattern='RData'), value=TRUE, ignore.case = TRUE)[1]
    print(sprintf("File used: %s", DATAFILENAME))
}else{
    stop(paste("No file with the following identifier '", opt$file, "'was found in \n", getwd()))
}
    
sTime = Sys.time()

SUBSET  = opt$subset

run_script = NULL
## Formating the sub-directory name
if(is.null(opt$output)){
    if(is.null(opt$subset)){
        subDir =  paste(toupper(opt$runtype),
            gsub('.RData', '', DATAFILENAME),
            format(sTime, "%d-%m-%Hh%M"), sep='-')
    }else{
        subDir =  paste(toupper(opt$runtype),
            opt$subset, gsub('.RData', '', DATAFILENAME),
            format(sTime, "%d-%m-%Hh%M"), sep='-')
    }
}else{
    subDir = opt$output
}
## adding sub-extension
subDir = paste0('./', subDir)
DATAFILENAME = paste0('../', DATAFILENAME)


## Choosing runscipt
if(grepl('uncertain', opt$runtype, ignore.case=TRUE)){
    run_script = '../main.10foldCV-uncertain.R'
    stop('this method is not yet implemented in this branch!')
}
if(grepl('(10fold|affinity|weighted|distonly)', opt$runtype, ignore.case =TRUE)){
    run_script = '../main.10foldCV.R'
    TYPE = toupper(opt$runtype)
}

if(grepl('ALL', opt$runtype, ignore.case = TRUE))
    run_script = '../main.all.R'

if(is.null(run_script) | is.null(subDir)){
    stop('something is wrong, the run script not found or no sub-directory is specified!')
}

## # Setting Hyper parameters
if(!is.null(opt$subset)){
    if(grepl('gmp', DATAFILENAME, ignore.case = TRUE))
        hyper = list(parasiteHyper =c(0.32, 1), hostHyper = c(0.94,1), etaSamplingSD = c(0.01)) 
    
    if(grepl('eid', DATAFILENAME, ignore.case = TRUE))
        hyper = list(parasiteHyper= c(0.35, 1), hostHyper =c(0.78, 1), etaSamplingSD = c(0.01)) 
    
}else{
    if(grepl('gmp', DATAFILENAME, ignore.case = TRUE)){
        hyper = list(parasiteHyper =c(0.25, 1), hostHyper = c(0.82,1), etaSamplingSD = c(0.005))
    }
    if(grepl('eid', DATAFILENAME, ignore.case = TRUE)){
        hyper = list(parasiteHyper =c(0.2, 1), hostHyper = c(0.73,1), etaSamplingSD = c(0.005))
    }
}


print("Arguments:")
print(opt)

source('library.R')
source('gen.R')
##Creating directory
reportFile <- "report.txt"
## Creating a subdirectory
dir.create(file.path(subDir))
## Setting the working directory
setwd(file.path(subDir))
## starting the process
sink(reportFile)
print(date())
print(DATAFILENAME)
print(run_script)
print(opt)
print(subDir)
## Process started at:
print(sprintf('Start time %s.',format(sTime, "%Y-%m-%d %H:%M:%S")))
##-----------------------------------------------------------
## Main running script
source(run_script, echo=TRUE, max.deparse.length=1e3)

##-----------------------------------------------------------
eTime = Sys.time()
print(sprintf('End time %s.',format(eTime, "%Y-%m-%d %H:%M:%S")))
## Processing time
print('Total processing time')
print(eTime - sTime)

######################
## Closing sink and reverting work directory.
sink()
system("grep '^[^>+;]' report.txt > report_clean.txt")
setwd('../')
if(!is.null(opt$email)){
    subj = '"End of sim "'
    body = paste(subDir)
    email = paste("echo '" ,body,"' | mailx -s ", subj, opt$email)
    system(email)
}
q('no')

