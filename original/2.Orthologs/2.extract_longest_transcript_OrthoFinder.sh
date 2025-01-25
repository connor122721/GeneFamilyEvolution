#!/usr/bin/env bash

# Load conda env 
module load anaconda/2023.07-py3.11
source activate msprime_env

# Loop through and extract the longest transcript in each protein fasta
for f in *faa; do python ~/OrthoFinder_source/tools/primary_transcript.py ${f}; done

# Rename and extract genes for European D. pulex
~/seqkit faidx \
--region-file /scratch/csm6hg/genomes/pulex_euro/primary_transcripts.txt  \
/scratch/csm6hg/genomes/primary/pulex.euro.protein.faa > \
/scratch/csm6hg/genomes/primary/pulex.euro.protein.sub.faa