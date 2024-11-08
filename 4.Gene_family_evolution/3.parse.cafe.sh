#!/usr/bin/env bash
# ijob -A berglandlab_standard -c 10 -p dev

# Working directory
cd /scratch/csm6hg/cafe_sep7/

# Extract significant gene families
echo $'#nexus\nbegin trees;'> run1_filt_hogs/Significant_trees.tre
grep "*" run1_filt_hogs/Base_asr.tre >> run1_filt_hogs/Significant_trees.tre
echo "end;">> run1_filt_hogs/Significant_trees.tre
awk '$2 < .05 {print $0}' run1_filt_hogs/Base_family_results.txt > run1_filt_hogs/Sig_at_p.05.txt

# Make HOG lists
grep "chitin" hogs.aug17.function.tsv > chitin.hog
grep "heat" hogs.aug17.function.tsv > heatshock.hog
grep "opsin" hogs.aug17.function.tsv > opsin.hog
grep "toll-" hogs.aug17.function.tsv > toll.hog
grep "immuno" hogs.aug17.function.tsv > immunoglob.hog
grep "vir" hogs.aug17.function.tsv > viral.hog
