#!/usr/bin/env bash
# ijob -A berglandlab -c 10 -p largemem 

# Load modules
module load gcc intel

# Cafe5 executable
cd ~
cafe=".local/bin/cafe5"

# Working directory
wd="/scratch/csm6hg/genomes/"

# Orthofinder output - HOGs
ortho=/scratch/csm6hg/cafe_sep7/hog_sep7_gene_counts.tsv

# Ultrametric time-calibrated tree from MCMCtree
tree="/scratch/csm6hg/cafe_aug28/mcmctree_busco_daphnia.nwk"

# Orthogroups from OrthoFinder
nohup ${cafe} \
-i ${ortho} \
-t ${tree} \
-o /scratch/csm6hg/cafe_sep7/run1_filt_hogs \
--cores 10 &

nohup ${cafe} \
-i ${ortho} \
-t ${tree} \
-o /scratch/csm6hg/cafe_sep7/run1_filt_k3_hogs \
-k 3 \
--cores 10 &

nohup ${cafe} \
-i ${ortho} \
-t ${tree} \
-o /scratch/csm6hg/cafe_sep7/run1_filt_k3_p_hogs \
-k 3 \
-p \
--cores 10 &

nohup ${cafe} \
-i ${ortho} \
-t ${tree} \
-o /scratch/csm6hg/cafe_sep7/run1_filt_e_hogs \
-e \
--cores 10 &

### Plotting ###

# Load conda env for speedseq
module load anaconda/2023.07-py3.11
source activate msprime_env
cd /scratch/csm6hg/cafe_aug18/

# plot using CAFE_fig
python3 ../scripts/CAFE_fig.py \
run1_filt_hogs/Base_report.cafe \
--dump plot/run1_filt_hogs -g .pdf 

cafeplotter \
-i run1_filt_hogs \
-o hogs_run1_plot \
--fig_width 5 \
--dpi 300 \
--expansion_color 'red' \
--contraction_color 'blue' \
--format 'pdf' \
--ignore_branch_length

cafeplotter \
-i run1_filt_k3_hogs \
-o hogs_run1_k3_plot \
--fig_width 5 \
--dpi 300 \
--expansion_color 'red' \
--contraction_color 'blue' \
--format 'pdf' \
--ignore_branch_length

cafeplotter \
-i run1_filt_k3_p_hogs \
-o hogs_run1_k3_p_plot \
--fig_width 5 \
--dpi 300 \
--expansion_color 'red' \
--contraction_color 'blue' \
--format 'pdf' \
--ignore_branch_length

