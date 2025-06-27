#!/usr/bin/bash -l
#SBATCH --job-name=snakemake_usearch
#SBATCH --output=job_%j_%x.out
#SBATCH --cpus-per-task=1
#SBATCH --time=0-06:00:00
#SBATCH --mem=1G

# Exit on first error and if any variables are unset
set -eu

conda activate snakemake_usearch

rm -rf logs tmp output
snakemake --unlock
snakemake --executor slurm --jobs 96 --notemp --profile profiles/biocloud
