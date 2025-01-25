#!/usr/bin/env bash
# ijob -A berglandlab_standard -c 10 -p largemem 

# Load conda env 
module load anaconda/2020.11-py3.8
source activate msprime_env
cd /scratch/csm6hg/genomes/primary/OrthoFinder/Results_Sep07/

# Annotate Orthogroups
annotate_orthogroups \
--orthogroups_tsv Phylogenetic_Hierarchical_Orthogroups/N0.tsv \
--hog True \
--fasta_dir /scratch/csm6hg/genomes/primary/ \
--file_endings faa \
--out hogs.sep7.function.tsv \
--simple Ture

# Make OrthoFinder Plots
orthofinder_plots \
--tree Species_Tree/SpeciesTree_rooted.txt \
--orthogroups_tsv Phylogenetic_Hierarchical_Orthogroups/N0.tsv \
--hog True \
--out orthofinder.sep7.plot
