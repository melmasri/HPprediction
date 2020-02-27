#!/usr/bin/env Rscript
library("optparse")

option_list = list(
	make_option(c("-m", "--model"), type="character",
                default="full", 
                help="possible models: full, distance, affinity [default= %default]."),
    make_option(c("-r", "--runtype"), type="character",
                default="CV", 
                help="possible run types: CV, all, uncertain [default= %default]."),
    make_option("--output", type="character",
                default=NULL, 
                help="output directory, the default is constructed from time stamp and other inputs."),
    make_option("--no.cycles", type="integer",
                default=20, 
                help="no. of cycles to run over the matrix [default= %default]."),
    make_option("--no.cores", type="integer",
                default=5, 
                help="number of cores to use (1, 2, ..) [default= %default]."),
    make_option("--data.file", type="character",
                default=NULL, 
                help="relative path to com RData file. File must contain a single HxJ interaction matrix named 'com' [default= %default]."),
    make_option("--alpha.rows", type="double",
                default=6, 
                help="hyperparameter for prior over rows affinity, effective under affinity and full models only [default= %default]."),
    make_option("--alpha.cols", type="double",
                default=0.03, 
                help="hyperparameter for prior over columns affinity, effective under affinity and full models only [default= %default]."),
    make_option("--cv.min.per.col", type="integer",
                default=1, 
                help="minimum number of interactions per column for cv, effective under runtype CV and uncertain only [default= %default]."),
    make_option("--seed", type="integer",
                default=23456, 
                help="integer value specifying starting seed [default= %default]."),
    make_option("--email", type="character",
                default=NULL,
                help="end-of-run email notification is sent to the supplied email if mailx is installed and configured")
); 

opt = parse_args(OptionParser(option_list=option_list))

##
set.seed(opt$seed)
## General variables
SAVE_PARAM = TRUE                    # should workspace be saved
SAVE_FILE = 'param.RData'            # name of output R workspace file
NO.CORES = opt$no.cores
SLICE = opt$no.cycles
PATH.TO.FILE = opt$data.file
ALPHA.ROWS = opt$alpha.rows
ALPHA.COLS = opt$alpha.cols
CV.MIN.PER.COL=opt$cv.min.per.col
MODEL = if(grepl('(full|aff|dist)', tolower(opt$model)))
            grep(tolower(opt$model), c('full', 'distance', 'affinity'), value=T) else NULL
COUNT = if(grepl('uncer', tolower(opt$runtype))) FALSE else TRUE
run_script = NULL


## validation test
if(grepl('dist', opt$model, ignore.case=TRUE) && grepl('(cv|uncer)', opt$runtype, ignore.case=TRUE) && CV.MIN.PER.COL==1){
    warning('minimum interactions per column are advised to be 2 when running cross-validation under the distance-only model. cv.min.per.col changed to 2!')
    CV.MIN.PER.COL = 2
}


## choosing an MCMC method 
if(grepl('CV', opt$runtype, ignore.case = TRUE))
    run_script = 'mainCV_network.R'

## if(grepl('latentnet', opt$runtype, ignore.case = TRUE))
##     run_script = 'extras/mainCV_network_latentnet.R'

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
    subDir =  paste(toupper(opt$runtype),toupper(opt$model),toupper(opt$output),
        format(sTime, "%d-%m-%Hh%M"), sep='-')
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
reportFile <- file(paste0(subDir, "report.txt"), 'wt')
    
## Setting the working directory

## starting the process
sink(reportFile)
sink(reportFile,type='message')
print(date())
print(run_script)
print(opt)
print(subDir)
## Git branch
print('Git info:')
system('git status -b -s',intern=TRUE)              # print branch info
system('git show --oneline -s',intern=TRUE)         # prin git commit 

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
sink(type='message')
close(reportFile)

system(paste("grep '^[^>+;]'", paste0(subDir, "report.txt"), ">", paste0(subDir, "report_clean.txt") ))
if(!is.null(opt$email)){
    subj = '"End of sim "'
    body = paste(subDir)
    email = paste("echo '" ,body,"' | mailx -s ", subj, opt$email)
    system(email)
}
q('no')

