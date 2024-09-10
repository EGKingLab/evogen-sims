#!/bin/bash
#SBATCH -p kingeg-lab,general,interactive,requeue
#SBATCH -A kingeg-lab
#SBATCH --job-name=spectral
#SBATCH --output=Combined_Spectral.out
#SBATCH --error=Combined_Spectral.err
#SBATCH --time=1-00:00:00 # Set the maximum runtime
#SBATCH -N 2
#SBATCH --cpus-per-task=128 # Number of CPU cores per task (adjust as needed)
#SBATCH --mem=490G # Total memory (adjust as needed)

module load r
module load quarto


quarto render /home/etb68/Mylab/evogen-sims/ET_Yeast/Writing1.dir/Code.dir/Combined_Spectral.qmd
