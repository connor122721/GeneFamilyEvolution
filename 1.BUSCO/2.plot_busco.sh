# #!/usr/bin/env bash

# Make plot for busco completeness across species
module load singularity

# Run 
singularity run /home/csm6hg/busco_v5.4.7_cv1.sif \
python3 /scratch/csm6hg/scripts/generate_plot.py \
-wd /scratch/csm6hg/genomes/proteins_species/primary_transcripts/summary