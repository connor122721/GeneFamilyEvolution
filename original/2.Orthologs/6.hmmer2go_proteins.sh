#!/usr/bin/env bash
#
#SBATCH -J hmmer2go # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH --mem 10G
#SBATCH -o /scratch/csm6hg/err/hmmer2go.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/err/hmmer2go.%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab

# Load modules
module load singularity
# singularity pull docker://sestaton/hmmer2go
# SLURM_ARRAY_TASK_ID=5

# WD
cd /scratch/csm6hg/genomes/proteins_species/
out="/scratch/csm6hg/genomes/proteins_species/hmmer_sep8"

# Protein fastas
pro=($(ls *faa))
samp=${pro[${SLURM_ARRAY_TASK_ID}]}
species=$( echo $samp | sed 's/.protein.faa//g' )
echo $species

# Run hmmer2go on each species' protein-fasta
singularity run /project/berglandlab/connor/hmmer2go_latest.sif \
hmmer2go run \
-i ${species}.protein.faa \
-d ~/PfamDB/Pfam-A.hmm \
-n 10 \
-o ${out}/${species}.genes_orf_Pfam-A.tblout

## Find Pfam -> GO term mapping
singularity run /project/berglandlab/connor/hmmer2go_latest.sif \
hmmer2go mapterms \
-i ${out}/${species}.genes_orf_Pfam-A.tblout \
-o ${out}/${species}.genes_orf_Pfam-A.map \
--map

## Create GAF file
singularity run /project/berglandlab/connor/hmmer2go_latest.sif \
hmmer2go map2gaf \
-i ${out}/${species}.genes_orf_Pfam-A.map \
-o ${out}/${species}.genes_orf_Pfam-A.map.gaf \
-s ${species}