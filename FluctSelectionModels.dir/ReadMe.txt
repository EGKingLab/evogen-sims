Code.dir ----> contains the SLiM codes for all selection models with:
	CS.dir for constant new optima selection, 
	NS.dir for null selection, 
	LinFS.dir for instantaneous selection, 
	SinFS.dir for gradual two equal seasons also known as simple gradual, 
	SinFSGen.dir for gradual four seasons, uneven season lengths and even distance to optima Δ, also known as complex gradual model I
	and lastly FourSeasFourAmp.dir for Gradual four unequal seasons, uneven season lengths and uneven distance to optima Δ also known as complex gradual model II. 
	
To run each model, make sure that you have created the directory for output files, and change directories as necessary.
You will need to run bash files interactive, but if necessary, you can create a slurm file to run the codes as well.

Output.dir -----> Contains folders for each model with log and csv files for allele frequency and phenotypes.

WritingFigures.dir -----> This folder is for the code used to produce figures in publication and their respective supplements.

The file "Scale_Estimation.qmd" shows how we estimated the phenotype scaling parameters

For some directories, you will find additional instructions in their corresponding readme files.

Note that all simulations fall under Wright-Fisher model.

