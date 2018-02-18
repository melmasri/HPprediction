#!/usr/bin/env Rscript
library("optparse")

option_list = list(
	make_option(c("-m", "--model"), type="character",
                default="full", 
                help="possible models: full, distance, affinity [default= %default]",
                metavar="character"),
    make_option(c("-r", "--runtype"), type="character",
                default="CV", 
                help="possible run types: CV, all, uncertain [default= %default]",
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
    make_option("--email", type="character",
                default=NULL,
                help="end-of-run email notification is sent to the supplied email if mailx is installed and configured",
                metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

##
set.seed(23456)
## General variables
SAVE_PARAM = TRUE                    # should workspace be saved
SAVE_FILE = 'param.RData'            # name of output R workspace file
NO.CORES = opt$no.cores
SLICE = opt$no.cycles
MODEL = if(grepl('(full|aff|dist)', tolower(opt$model)))
            grep(tolower(opt$model), c('full', 'distance', 'affinity'), value=T) else NULL
COUNT = if(grepl('uncer', tolower(opt$runtype))) FALSE else TRUE
run_script = NULL

## choosing an MCMC method 
if(grepl('CV', opt$runtype, ignore.case = TRUE))
    run_script = 'mainCV_network.R'

if(grepl('latentnet', opt$runtype, ignore.case = TRUE))
    run_script = 'mainCV_network_latentnet.R'

if(grepl('uncer', opt$runtype, ignore.case = TRUE))
    run_script = 'mainCV-uncertain_network.R'

if(grepl('all', opt$runtype, ignore.case = TRUE))
    run_script = 'main_network.R'

## Setting up the file name 
sTime = Sys.time()
## Formating the sub-directory name
if(is.null(opt$output)){
    subDir =  paste(toupper(opt$runtype),toupper(opt$model), 
        format(sTime, "%d-%m-%Hh%M"), sep='-')
}else{
    subDir = opt$output
}
## adding sub-extension
subDir = paste0('./', subDir, '/')

if(is.null(MODEL) | is.null(run_script)){
    stop('something is wrong, run script or model not specified!')
}


print("Arguments:")
print(opt)

source('network_analysis.R')
source('network_MCMC.R')
## Creating a subdirectory
dir.create(file.path(subDir))
## report file
reportFile <- paste0(subDir, "report.txt")
## Setting the working directory

## starting the process
sink(reportFile)
print(date())
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
system(paste("grep '^[^>+;]'", reportFile, ">", paste0(subDir, "report_clean.txt") ))

if(!is.null(opt$email)){
    subj = '"End of sim "'
    body = paste(subDir)
    email = paste("echo '" ,body,"' | mailx -s ", subj, opt$email)
    system(email)
}
q('no')

