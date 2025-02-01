#!/usr/bin/env bash
#
#SBATCH -J nextflow_mega # Job name
#SBATCH --ntasks-per-node=4 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # days
#SBATCH --mem 50G
#SBATCH -o /project/berglandlab/connor/err/nextflow.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/nextflow.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Modules to load
module load nextflow

# Run nextflow
# nextflow run extract_busco_genes.nf -profile slurm -resume
nextflow run main.nf -profile slurm -resume
