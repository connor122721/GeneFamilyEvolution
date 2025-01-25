#!/usr/bin/env bash
#
#SBATCH -J Busco # A single job name for the array
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00 # 1 hour
#SBATCH --mem 30G
#SBATCH -o /scratch/csm6hg/err/busco.new.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/err/busco.new.%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab_standard

# Docker pull busco - run once
module load singularity
# singularity pull docker://ezlabgva/busco:v5.4.7_cv1
# SLURM_ARRAY_TASK_ID=6

# Move to data directory
cd /scratch/csm6hg/genomes/proteins_species/primary_transcripts/

# Protein fastas
pro=($(ls *faa))
samp=${pro[${SLURM_ARRAY_TASK_ID}]}
species=$( echo $samp | sed 's/.protein.faa//g' )
echo $species

# Run Busco
singularity run /home/csm6hg/busco_v5.4.7_cv1.sif \
busco \
-i /scratch/csm6hg/genomes/proteins_species/primary_transcripts/${samp} \
-c 5 \
--out_path /scratch/csm6hg/genomes/proteins_species/primary_transcripts/ \
-l arthropoda_odb10 \
-o ${species} \
-m proteins
