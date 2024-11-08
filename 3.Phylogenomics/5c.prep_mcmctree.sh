#!/usr/bin/env bash

# Make prep tree for MCMCtree
module load anaconda/2023.07-py3.11
source activate msprime_env

cd /scratch/csm6hg/genomes/sco_aug28/trees

# Add time constraints
python /scratch/csm6hg/scripts/mcmctree_prep.py \
--left_species pulex.nam \
--right_species magna \
--lower_bound 130 \
--upper_bound 150 \
--tree sco_species_astralpro \
| python /scratch/csm6hg/scripts/mcmctree_prep.py \
--left_species magna \
--right_species sinensis \
--lower_bound 21.5 \
--upper_bound 22.4 \
--tree - \
--add_header \
> sco_daphnia_mcmc_1.19.nwk

# Move to MCMCtree folder
# mkdir /scratch/csm6hg/genomes/sco_aug28/mcmctree_aug28
cp /scratch/csm6hg/genomes/sco/trees/sco_daphnia_mcmc.nwk \
/scratch/csm6hg/genomes/sco_aug28/mcmctree_aug28

# Make mega sequence alignment for MCMCtree
# Extarct each SCO - FOR MCMCTree
while read -r line;
    do
      # line=101491at6656; spp='magna'
      echo ${line}
            
      # Rename
      sed "s/|${line}//g" /scratch/csm6hg/genomes/sco_aug28/align/${line}.align.fa.clipkit > \
      /scratch/csm6hg/genomes/sco_aug28/spp/${line}.align.clip.fa

done < /scratch/csm6hg/genomes/sco_aug28/sco.aug28.txt

# Concatenate aligned genes
~/seqkit concat /scratch/csm6hg/genomes/sco_aug28/spp/*.fa  \
-o /scratch/csm6hg/genomes/sco_aug28/spp/sco_daphnia_mega_clip.fa

cp /scratch/csm6hg/genomes/sco_aug28/spp/sco_daphnia_mega_clip.fa \
/scratch/csm6hg/genomes/sco_aug28/mcmctree_aug28/sco_daphnia_mega_clip
