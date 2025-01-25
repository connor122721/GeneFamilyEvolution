#!/usr/bin/env bash
#
#SBATCH -J TREE # A single job name for the array
#SBATCH --ntasks-per-node=4 # one core
#SBATCH -N 1 # on one node
#SBATCH --mem 50G
#SBATCH -o /scratch/csm6hg/err/maketree.out # Standard output
#SBATCH -e /scratch/csm6hg/err/maketree.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Run PAML package MCMCtree
export PATH=~/paml-4.10.7/bin/:$PATH

# Working directory
cd /scratch/csm6hg/genomes/sco_aug28/mcmctree_aug28_run3

# Run mcmctree
mcmctree mcmctree.ctl
