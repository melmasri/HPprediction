### Global settings
rm(list=ls())                           # clear environment
devtools::load_all()

set.seed(23456)                         # fixing seed for replication results

## General variables
SAVE_PARAM = TRUE                    # should workspace be saved
SAVE_FILE = 'param.RData'            # name of output R workspace file
NO.CORES = 1                         # no. of cores effective for runtype CV and uncertain
SLICE = 200                           # no. of cycles over the matrix
PATH.TO.FILE = NULL                  # relative path to com RData file., NULL to execute '   source('example-GMPD/load_GMPD.R')'
ALPHA.ROWS = 6                          # hyperparameter for prior over rows (>0), effective under affinity and full models
ALPHA.COLS = 0.03                       # hyperparameter for prior over cols (>0), effective under affinity and full models
CV.MIN.PER.COL=1                     # minimum minimum number of interactions per column for cv, effective under runtype CV and uncertain only.
MODEL = 'full'                       # model type options {'full', 'distance', 'affinity'}
COUNT = FALSE                         # only needed to be FALSE for uncertainty subsetting by years

subDir = ''                             #sub-directory where the saved output and other results are desired, directory must be created

### creating a directory if it doesn't exists
if(subDir!='')  dir.create(file.path(subDir))

## loading mammal supertree included in Fritz et al. (2009)
source('data/download_tree.R')       # see variable 'tree'

## loading GMPD
if(exists("PATH.TO.FILE") && !is.null(PATH.TO.FILE)){
    if(grepl('.rds', PATH.TO.FILE, ignore.case = TRUE))
        com <- readRDS(PATH.TO.FILE) else 
    load(PATH.TO.FILE)
}else{
    source('data/load_GMPD.R')           # see matrix 'com'    
}

#################################
## Running the possible 3 options 

## 1) Running the model on all the data without CV
source('main_network.R')

## 2) Running CV script
source('mainCV_network.R')

## 3) Running uncertainty script with CV
source('mainCV-uncertain_network.R')


