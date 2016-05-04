## what to run
## Rscript runscript.R runtype datatype subtype

## runtype: uncertain, 10fold, ALL, nodist, weighted, NN.
## subtype: subset, anything
## datatype: GMP EID


# A general running script
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
##library scripts
if(length(args)<3){
    print('len args min 3.')
    q('no')
}

runtype = args[1]
datatype = args[2]
subtype  = args[3]
sTime = Sys.time()

if(grepl('GMP',datatype)){
    DATAFILENAME = if(grepl('subset', subtype)| grepl('uncertain', runtype))'../comGMPD-year.RData' else '../comGMPD.RData'
}
if(grepl('EID',datatype)){
    DATAFILENAME = if(grepl('subset', subtype)| grepl('uncertain', runtype)) '../comEID-subset.RData' else  '../comEID-PS.RData'
}

SUBSET = if(grepl('subset', subtype)) TRUE else FALSE

if(grepl('GMP', DATAFILENAME)){
    dataset='gmp'
    subset= if(SUBSET) 'Carnivora-' else ''
}
if(grepl('EID', DATAFILENAME)){
    dataset='eid'
    subset= if(SUBSET) 'Rodentia-' else ''
}

run_script = NULL
subDir = NULL
if(grepl('uncertain', runtype)){
    run_script = '../main.10foldCV-uncertain.R'
    subDir <- paste0("./Uncertain-10foldCV-",subset, dataset, format(sTime, "-%d-%m-%Hh%M"))
}

if(grepl('10fold', runtype)){
    run_script = '../main.10foldCV.R'
    subDir <- paste0("./10foldCV-",subset, dataset, format(sTime, "-%d-%m-%Hh%M"))
}

if(grepl('ALL', runtype)){
    run_script = '../main.all.R'
    subDir <- paste0("./",subset, dataset, format(sTime, "-%d-%m-%Hh%M"))
}

if(grepl('nodist', runtype)){
    run_script = '../main.10foldCV-nodist.R'
    subDir <- paste0("./nodist-10foldCV-",subset, dataset, format(sTime, "-%d-%m-%Hh%M"))
}

if(grepl('weighted', runtype)){
    run_script = '../main.10foldCV-weighted.R'
    subDir <- paste0("./Weighted-10foldCV-",subset, dataset, format(sTime, "-%d-%m-%Hh%M"))
}

if(grepl('NN', runtype)){
    run_script = '../main.NN.R'
    subDir <- paste0("./NN-10foldCV-",subset, dataset, format(sTime, "-%d-%m-%Hh%M"))
}

if(is.null(run_script) | is.null(subDir)){
    print('something is wrong')
    q('no')
}

title = paste0(dataset,' ',subset,' all data')
source('library.R')
##source('gen.yEta.R')
source('gen.R')
## source('gen.wEta.R')
##Creating directory

reportFile <- "report.txt"
## Creating a subdirectory
dir.create(file.path(subDir))
## Setting the working directory
setwd(file.path(subDir))
## starting the process
sink(reportFile)
print(title)
print(date())
print(DATAFILENAME)
print(run_script)
print(args)
## Process started at:
print(sprintf('Start time %s.',format(sTime, "%Y-%m-%d %H:%M:%S")))
##-----------------------------------------------------------
## print(sprintf('lambda is %0.3f and eta is %0.3f', lambda_phy, eta))
## Main running script
source(run_script, echo=TRUE, max.deparse.length=1e3)
##source('../main.100sim-10foldCV.R', echo=TRUE)
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
subj = '"End of sim "'
body = paste(title, subDir)
email = paste("echo '" ,body,"' | mailx -s ", subj, " elmasri.m@gmail.com")
system(email)
q('no')

