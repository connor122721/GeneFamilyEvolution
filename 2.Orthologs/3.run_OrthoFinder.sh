#!/usr/bin/env bash
#
#SBATCH -J orthofinderCrust # A single job name for the array
#SBATCH --ntasks-per-node=15 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 day
#SBATCH --mem 50G
#SBATCH -o /scratch/csm6hg/err/ortho.crust.new.out # Standard output
#SBATCH -e /scratch/csm6hg/err/ortho.crust.new.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load conda env 
module load anaconda/2023.07-py3.11
source activate msprime_env

# Run OrthoFinder on primary protein transcripts
cd /scratch/csm6hg/crust.genomes/primary/
~/OrthoFinder_source/orthofinder.py -f . -t 4