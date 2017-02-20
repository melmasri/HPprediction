The `ICM` branch fits the model using the "Iterative Conditional Modelling" approach, discussed Besag(74) Section 6.1, which is called the "Coding" method.


# General setting for script running
The `R` script called `runscript.R` fits the model on some data input. In the terminal `runscript.R` can be executed with some arguments with the command `Rscript`.
The input arguments are:
```
Usage: runscript.R [options]


Options:
	-f CHARACTER, --file=CHARACTER
		dataset file name, or an identifiable portion of it, if multiple files the first one is used. The file should have the RData extension.

	-r CHARACTER, --runtype=CHARACTER
		possible run types, ALL, 10fold, affinity, distonly, weighted [default= 10fold]

	-s CHARACTER, --subset=CHARACTER
		possible subset choices, for GMP it is Carnivora, for and EID it is Rodentia.

	--output=CHARACTER
		output directory, the default is constructed from time stamp and other inputs.

	--no.cycles=INTEGER
		no. of cycles to run over the matrix [default= 20].

	--no.cores=INTEGER
		number of cores to use (1, 2, ..) [default= 5].

	--icm.horiz=LOGICAL
		block updates are done horizontally (ICM.HORIZ=TRUE), else diagonally [default= TRUE].

	--email=CHARACTER
		end-of-run email notification is sent to the supplied email if mailx is installed and configured

	-h, --help
		Show this help message and exit
```

+ `runtype`: 
    + `ALL`: to fit the model on the whole dataset without CV.
    + `10fold`: to fit the full model using 5-fold CV.
    + `affinity`: to fit the affinity-only model using 5-fold CV.
    + `distonly`: to fit the dist-only model using 5-fold CV.
    + `weighted`: to fit the full model using 5-fold CV and the weighted (count) interaction matrix.
+ `no.cycles`: two methods to calculate the number of MCMC samples given the `no.cycles`
    + with diagonal block update `icm.horiz=FALSE`, the model cycles over the number of columns in the interaction matrix. For a 50x100 matrix, a cycle has 50 samples for each column parameter, and 100 for each row parameter.
    + with horizontal block update `icm.horiz=TRUE`, the model cycles over the number of rows in the interaction matrix. For a 50x100 matrix, a cycle has 50 samples for each column parameter, and 1 for each row parameter.


For example to run the 10fold affinity only model using a data set file called `Dataset123.RData`,

`
Rscript runscript.R -f Dataset123 -r affinity
`

A an output folder is created with the following name `AFFINITY-Dataset123-TIMESTAMP`.
      

    
