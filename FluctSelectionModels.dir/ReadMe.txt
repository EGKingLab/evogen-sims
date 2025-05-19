Code.dir ----> contains the SLiM codes for all selection models with:
	CS.dir for constant selection, 
	NS.dir for neutral selection, 
	LinFS.dir for instantaneous selection, 
	SinFS.dir for gradual two equal seasons, 
	and lastly SinFSGen.dir for Gradual four unequal seasons selection models. 

To run each model, make sure that you have created the directory for output files, and change directories as necessary.
You will need to run bash files interactive, but if necessary, you can create a slurm file to run the codes as well.

Note that all simulations fall under Wright-Fisher model

Output.dir -----> Contains folders for each model with log and csv files for allele frequency and phenotypes.

WritingFigures.dir -----> This folder is for the code used to produce figures in publication and their respective supplements.

The file "Scale_Estimation.qmd" shows how we estimated the phenotype scaling parameters
