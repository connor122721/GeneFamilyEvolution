#!/usr/bin/env bash
#
#SBATCH -J sco_spp_tree # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 days
#SBATCH --mem 2G
#SBATCH -o /scratch/csm6hg/err/sco.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/err/sco.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# Make full gene tree
cd /scratch/csm6hg/genomes/sco_aug28/trees
cat *.tre.treefile > busco_sco_genes.tre

# Get rid of sco geneid
while read -r line;
do
    # line=101491at6656
    echo $line
    sed -i 's/'${line}'//g' busco_sco_genes.tre
done < /scratch/csm6hg/genomes/sco_aug28/sco.aug28.txt

# Remove "|" character
sed -i 's/|:/:/g' busco_sco_genes.tre
sed -i 's/|2:/:/g' busco_sco_genes.tre

# Run Astral-pro to create species tree
~/ASTER-1.15/bin/astral-pro \
-i busco_sco_genes.tre \
-R \
--root "carinata" \
-o sco_species_astralpro > \
astralpro.log

# Open density tree with beast2
/home/csm6hg/beast/bin/densitree

# Make prep tree for MCMCtree
module load anaconda/2020.11-py3.8
source activate msprime_env
python /scratch/csm6hg/scripts/mcmctree_prep.py \

python mcmctree_tree_prep.py \
--left_species carinata \
--right_species magna \
--lower_bound 104 \
--upper_bound 1 \
--tree input.nwk \
| python mcmctree_tree_prep.py \
--left_species Aquilegia_coerulea \
--right_species Arabidopsis_thaliana \
--lower_bound 130 \
--upper_bound 130 \
--tree - \
--add_header \
> output.nwk