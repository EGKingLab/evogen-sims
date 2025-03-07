#!/bin/bash
#SBATCH -p general,interactive,requeue
#SBATCH -A general
#SBATCH --job-name=NS Slim
#SBATCH --output=NS.out
#SBATCH --error=NS.err
#SBATCH --time=7:00:00 # Set the maximum runtime
#SBATCH --cpus-per-task=56 # Number of CPU cores per task (adjust as needed)
#SBATCH --mem=56G # Total memory (adjust as needed)

module load slim

bash BashNS_Par4_Nivalis.sh
