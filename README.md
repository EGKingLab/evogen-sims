# evogen-sims

This is a repository for "Dissecting fluctuating selection: A unified population and quantitative genetics" project.

The major folders are:
1. FluctSelectionModels.dir - For models of selection and constant population size
2. FluctPopulationModels.dir - For models of fluctuating population size and selection (Optima)

Each directory and its child directories are described below:

├── evogen_sims.Rproj

├── FluctSelectionModels.dir: This folder contains all simulations of different types of selection and constant population size

│   ├── Code.dir: This subdirectory with child directories simulates different scenarios of selection and constant population size as described below 

- Note that each child folder has a bash file with required parameters and their combinations, and it takes a slim file to run simulations.
- The slim file takes each parameter combination, runs the evolutionary simulation and writes the output files to the corresponding folder.
- Therefore, the user must first ensure that the output folders are created.
- Furthermore, they should make sure that they have requested enough resources before running some of these simulations, especially those that require extended generations.
- Although we run our simulations using bash files, one can also use Slurm

│   ├── Code.dir

│   │   ├── NS.dir: Null selection models with constant population size

│   │   │ ├── BashNS_Mod_RandSeed.s -->The bash file with parameters

│   │   │ ├── NS_Mod.slim --> The slim file to run simulations

│   │   ├── CS.dir: Constant selection and population size

│   │   │ ├── BashCS_Mod_RandSeed.sh --> The bash file with parameters

│   │   │ ├── CS_Mod.slim --> The slim file to run simulations

│   │   ├── LinFS.dir: Instantaneous optima change with two equal seasons

│   │   │ ├── BashFS_Mod_RandSeed.sh --> The bash file with parameters

│   │   │ ├── FS_Mod.slim --> The slim file to run simulations

│   │   ├── SinFS.dir: Gradual optima change with two equal seasons

│   │   │ ├── BashGradI_Mod_RandSeed.sh --> The bash file with parameters

│   │   │ ├── GradI_Mod.slim --> The slim file to run simulations

│   │   └── SinFSGen.dir: Gradual optima change with even optima and uneven season length

│   │   │ ├── Bash_UpGradII_Mod_RandSeed.sh --> The bash file with parameters

│   │   │ ├── UpdatedGradII_Mod.slim --> The slim file to run simulations

│   │   ├── FourSeasFourAmp.dir: Gradual optima change with uneven optima and uneven season length

│   │   │ ├── BashFourSeasFourAmp.sh --> The bash file with parameters

│   │   │ ├── FourSeasFourAmp.slim --> The slim file to run simulations

│   │   │ ├──

│   │   ├── ExtendedGenMuRatePopLinFS.dir: Extended models (Generation = 10N) with mutation rate and constant population size

│   │   │ ├── Bash_ExtendedGenMuRatePopLinFS.sh --> The bash file with parameters

│   │   │ ├── ExtendedGenMuRatePopLinFS.slim --> The slim file to run simulations

│   │   │ ├──

│   │   ├── Scale_Estimation.qmd --> This is the code we used to estimate the constant C that influences the steepness of the fitness surface

│   ├── Output.dir: The main output directory with child directories. Make sure to create them before running the simulations.

│   │   ├── CS.dir

│   │   ├── ExtendedGenMuRatePopLinFS.dir

│   │   ├── FourSeasFourAmp.dir

│   │   ├── LinFS.dir

│   │   ├── LongRuns.dir

│   │   ├── NS.dir

│   │   ├── SinFS.dir

│   │   └── SinFSGen.dir

│   ├── ReadMe.txt--> ReadMe file for fluctuating selection only

│   └── WritingFigures.dir: The following folder contains the code and the directories for the expected output of the published figures and supplements.

│       ├── CodeFigures.dir - The code folder. When running the code, make sure you have loaded all required packages, including Quarto and requested enough resources.
                              The number of cores requested for parallelization can be adjusted based on the institution's resources. 

│       ├── OutputFigures.dir - The output folder for figures.

│       └── ReadMe.txt --> ReadMe file for data processing/ Figures production only

├── FluctPopulationModels.dir: Different models selection and fluctuating population size

│   │   ├── ConstSelInstatPop: A folder for fluctuating population sizes and constant selection models.

│   │   │ ├── BashConstSelInstatPop.sh - a bash file with considered parameters and their combination

│   │   │ ├──ConstSelInstatPop.slim - A slim file that takes in the parameters from the bash file and writes the output files in the designated folder

│   │   ├── InstSelInstatPop: A folder for an instantaneous selection and population size

│   │   │ ├── BashInstSelInstatPop.sh - The bash file with parameters

│   │   │ ├── InstSelInstatPop.slim - The slim file to run simulations

│   │   └── NeutrSelInstPop: Neutral selection with instantaneous change in population size

│   │   │ ├── BashNeutrSelInstatPop.sh - The bash file with parameters

│   │   │ ├──  NeutrSelInstatPop.slim - The slim file to run simulations

│   ├── Output.dir: The corresponding output folders should be created before running the simulations for each scenario

│   │   ├── ConstSelInstatPop

│   │   ├── InstSelInstatPop

│   │   └── NeutrSelInstatPop

│   ├── ProjectNotes.qmd

│   └── ReadMe.txt

└────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

