#!/usr/bin/env bash
#
#SBATCH -J hyphy_mega # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-02:00 # 2 hours
#SBATCH --mem 2G
#SBATCH -o /scratch/csm6hg/err/hyphy.gene.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/err/hyphy.gene.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# Working directory: SLURM_ARRAY_TASK_ID=1
cd /scratch/csm6hg/codeml/

# Where protein SCO MSAs are: 19308
sco="/scratch/csm6hg/genomes/primary/OrthoFinder/Results_Sep07/Orthogroup_Sequences"

# Conditional make mega list for codeml
if [ -f mega.list.OG ]; then
  echo "Mega list is there."
else
  echo "Empty mega list."
  ls -d ${sco}/*fa > \
  mega.list.OG
fi

# The gene we are working on; grep "OG0000132" -n mega.list.OG
gene=$( sed -n "${SLURM_ARRAY_TASK_ID}p" mega.list.OG2 )
id=$( sed -n "${SLURM_ARRAY_TASK_ID}p" mega.list.OG2 | rev | cut -f1 -d "/" | rev | sed 's/.fa//')
echo $gene "Orthogroup =" $id

# Genomes
genome="/scratch/csm6hg/genomes"

# Conditional make mega list for codeml
if [ -d mega_OG ]; then
  echo "WD there."
else
  echo "Empty mega WD."
  mkdir mega_OG
fi

# Conditional make gene working directory
if [ -d mega_OG/${id} ]; then
  echo "WD there."
else
  echo "Empty mega WD."
  mkdir mega_OG/${id}
fi

# Move into working directory
cd mega_OG/${id}

#######################################################
################ Arrange codon MSAs ###################
#######################################################

# Make list of genes
cat ${gene} | grep ">" | sed 's/>//g' > ${id}.genes

# Extract gene sequences
for i in ${genome}/*/*cds_from_genomic_rename.fna; do
   echo $i
   
   # Extract gene sequences (DNA)
   ~/seqkit faidx ${i} -l ${id}.genes >> ${id}.fa

done

# Align the nuc
module load gcc/11.4.0 openmpi/4.1.4 mafft/7.505
mafft --auto \
${id}.fa > \
${id}.align.fa

# Align aa
mafft --auto \
${gene} > \
${id}.align.faa

# Download pal2nal: http://www.bork.embl.de/pal2nal/#RunP2N
# Translate and make input for codeml
~/pal2nal.pl \
${id}.align.faa \
${id}.align.fa \
-output fasta -nogap -nomismatch > \
${id}.trans.fa

# Clip out poorly aligned regions
module load anaconda/2023.07-py3.11
source activate papaml
module load apptainer
apptainer run ~/gblocks_latest.sif Gblocks \
${id}.trans.fa \
-t=c \
-p=n \
-b5=h \
-e=.gb

# Make fasta > 1 line per sequence
perl ../../one_line_fasta.pl ${id}.trans.fa.gb

# Remove spaces introduced through Gblocks processing
cat ${id}_one_line.fa | tr -d " " > ${id}.trans.clip.fa

# Run iqtree2 for a gene tree
~/iqtree2 \
-redo \
-s ${id}.trans.clip.fa \
-T 1 \
--prefix ${id}.tre

# Remove dup sequences
hyphy ~/hyphy-analyses/remove-duplicates/remove-duplicates.bf \
--msa ${id}.trans.clip.fa \
--tree "/scratch/csm6hg/codeml/mega_OG/${id}/${id}.tre.treefile" \
--output ${id}.uniques.nxh

# Conditional run Hyphy
# Run aBSREL: install hyphy -> conda install -c bioconda hyphy
if [ -f ${id}.uniques.nxh ]; then
  echo "Duplicates found, running hyphy."
  hyphy absrel \
  --alignment ${id}.uniques.nxh \
  --tree ${id}.tre.treefile \
  --output ${id}.ABSREL.json
else
  echo "No duplicates found."
  hyphy absrel \
  --alignment ${id}.trans.clip.fa \
  --tree ${id}.tre.treefile \
  --output ${id}.ABSREL.json
fi

# Remove intermediate files
rm ${id}*_one_line.fa
rm ${id}*.splits.nex
rm ${id}*.model.gz
rm ${id}*.varsites.phy
rm ${id}*.uniqueseq.phy
rm ${id}*.iqtree
rm ${id}*.bionj
rm ${id}*.contree
rm ${id}*.mldist
rm ${id}*.ufboot
rm ${id}*.ckp.gz
rm ${id}*.log

# Finish gene
echo "Finish gene: " $id
date