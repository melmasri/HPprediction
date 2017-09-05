Automation of the implemented model through `Rscript` command line.

# General setting for  running

The `R` script called `runscript.R` fits the 3 model variations (full, distance, affinity) using either cross-validation (CV) or the full dataset. In the terminal `runscript.R` can be executed with some arguments with the command `Rscript` with the help of `optparse` library. The input arguments are:
```
Usage: extras/runscript.R [options]


Options:
	-m CHARACTER, --model=CHARACTER
		possible models: full, distance, affinity, uncertain, [default= full]

	-r CHARACTER, --runtype=CHARACTER
		possible run types: all, CV, [default= CV]

	--output=CHARACTER
		output directory, the default is constructed from time stamp and other inputs.

	--no.cycles=INTEGER
		no. of cycles to run over the matrix [default= 20].

	--no.cores=INTEGER
		number of cores to use (1, 2, ..) [default= 5].

	--email=CHARACTER
		end-of-run email notification is sent to the supplied email if mailx is installed and configured

	-h, --help
		Show this help message and exit
```

+ `runtype`: 
    + `all`: to fit the model on the whole dataset without CV, see `main_network.R` for more details;
    + `CV` : to fit the model using cross-validation, requires library `parallel`, see `mainCV_network.R` for more details.
+ `model`
    + `full`: to fit the full model;
    + `distance`: to fit the distance model;
    + `affinity`: to fit the affinity model.
+ `no.cycles`: the number of MCMC iteration to run, including burn in:
    + for the affinity model, since the full joint posterior is used, this corresponds to the number of MCMC samples;
    + for full and distance models, since the conditional row joint posterior is used, this corresponds to the number of cycles over the rows of the interaction matrix. An interaction matrix of `H` rows requires `H` cycles for a single MCMC sample. 
+ `no.cores`: the number of cores to use in the `CV` runtype.
+ `email`: a notification sent to the specified email is a SMTP server is setup on the local machine.  


# Example
To run the full model with cross validation: navigate the main directory (`HP-prediction`), and run the following command on the shell. 

`
Rscript extras/runscript.R -m full -r cv --no.cycles= 1000 --no.cores 2
`

A an output folder will be created with the following name `CV-FULL-TIMESTAMP`. 
      
# Note
The automation procedure uses the scripts `main_network.R` and `mainCV_network.R` in the main directory. Thus one will need to comment out the general setting variables in each of those two scripts. Those variables are:

```R
## General variables
## please specify the following parameters
SAVE_PARAM = TRUE                    # should workspace be saved
SAVE_FILE = 'param.RData'            # name of output R workspace file
MODEL = 'full'                       # full, distance or affinity
SLICE = 100                          # no of iterations
subDir = ''                          # directory to print the results 
NO.CORES = 2                         # maximum cores to use
```
