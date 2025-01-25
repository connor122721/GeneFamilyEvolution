#!/usr/bin/env bash
#
#SBATCH -J sco_trees # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --mem 2G
#SBATCH -o /scratch/csm6hg/err/trees/sco.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/err/trees/sco.%A_%a.err # Standard error
#SBATCH -p dev
#SBATCH --account berglandlab_standard

# Generate tree for SCO genes
# SLURM_ARRAY_TASK_ID=10
line=$( sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch/csm6hg/genomes/sco_aug28/sco.aug28.txt )
echo $line

# Run iqtree2
~/iqtree2 \
-redo \
-bb 1000 \
-s /scratch/csm6hg/genomes/sco_aug28/align/$line.realign.clip.fa \
-T 1 \
--prefix /scratch/csm6hg/genomes/sco_aug28/trees/$line.tre

# Remove temp files
rm /scratch/csm6hg/genomes/sco_aug28/trees/$line.tre*splits.nex
rm /scratch/csm6hg/genomes/sco_aug28/trees/$line.tre*model.gz
rm /scratch/csm6hg/genomes/sco_aug28/trees/$line.tre*varsites.phy
rm /scratch/csm6hg/genomes/sco_aug28/trees/$line.tre*uniqueseq.phy
rm /scratch/csm6hg/genomes/sco_aug28/trees/$line.tre*iqtree
rm /scratch/csm6hg/genomes/sco_aug28/trees/$line.tre*bionj
rm /scratch/csm6hg/genomes/sco_aug28/trees/$line.tre*contree
rm /scratch/csm6hg/genomes/sco_aug28/trees/$line.tre*mldist
rm /scratch/csm6hg/genomes/sco_aug28/trees/$line.tre*ufboot
rm /scratch/csm6hg/genomes/sco_aug28/trees/$line.tre*ckp.gz
rm /scratch/csm6hg/genomes/sco_aug28/trees/$line.tre*log

# Finish
echo "Finish tree:" ${line}
