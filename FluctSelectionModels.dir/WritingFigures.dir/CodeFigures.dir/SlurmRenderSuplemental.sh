#!/bin/bash
#SBATCH -p requeue,logical_cpu
#SBATCH -A kingeg-lab
#SBATCH --job-name=quarto_suppfigs
#SBATCH --output=quarto_suppfigs_%j.out
#SBATCH --error=quarto_suppfigs_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=490G

set -eo pipefail

CONDA_BASE="/cluster/software/SPACK/SPACK_v0.20_dev_a2/spack/opt/spack/linux-almalinux8-x86_64/gcc-12.3.0/miniconda3-4.10.3-c6moxpqnii2vbelazwz5onnnnsh3cbzm"

source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate r-env

export QUARTO_R="$CONDA_PREFIX/bin/R"

export PATH="$CONDA_PREFIX/bin:$PATH"

export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"

echo "Job started on $(hostname)"

echo "Date: $(date)"

echo "CONDA_BASE=$CONDA_BASE"

echo "CONDA_PREFIX=$CONDA_PREFIX"

echo "R:"

which R

R --version

echo "Rscript:"

which Rscript

echo "Quarto:"

which quarto

quarto --version

echo "R libraries:"

Rscript -e '.libPaths(); cat("rlang path:", find.package("rlang"), "\n"); cat("rlang version:", as.character(packageVersion("rlang")), "\n"); cat("tidyr path:", find.package("tidyr"), "\n"); cat("tidyr version:", as.character(packageVersion("tidyr")), "\n")'

quarto render SupplementalFigures.qmd --to docx --verbose --no-cache
