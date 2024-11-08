# #!/usr/bin/env bash

# Working directory
cd /scratch/csm6hg/genomes/sco
module load gcc/9.2.0 openmpi/3.1.6 mafft/7.475

for f in /scratch/csm6hg/genomes/busco/*.sco.txt
do

    # f="/scratch/csm6hg/genomes/busco/pulex.euro.sco.txt"

    echo ${f}
    spp=$( echo ${f} | cut -f6 -d "/" | sed 's/.sco.txt//g' )
    echo ${spp}

    # Extract SCO from genome file
    ~/seqkit faidx \
    -l ${f} \
    /scratch/csm6hg/genomes/proteins_species/${spp}*.faa > \
    /scratch/csm6hg/genomes/sco_aug28/${spp}.sco.faa

done

# After this step run: rename SCO fasta script (2.extract_busco_genes.R)

# Go through each spp
for f in /scratch/csm6hg/genomes/sco_aug28/rename/*sco.faa 
do
    echo ${f}
    spp=$( echo ${f} | cut -f7 -d "/" | sed 's/.rename.sco.faa//g' )
    echo ${spp}

        # Extarct each SCO 
        while read -r line;
        do
            # line=101491at6656; spp='magna'
            echo ${line}
            
            # Extract sequence
            ~/seqkit faidx ${f} ${line} > /scratch/csm6hg/genomes/sco_aug28/genes/${line}.${spp}.fa

            # Rename
            sed -i "s/>${line}/>${spp}|${line}/g" /scratch/csm6hg/genomes/sco_aug28/genes/${line}.${spp}.fa
        done < /scratch/csm6hg/genomes/sco_aug28/sco.aug28.txt

done

# Align each SCO with mafft
while read -r line;
do
    # line=101491at6656
    echo ${line}

    cat /scratch/csm6hg/genomes/sco_aug28/genes/*${line}*.fa > \
    /scratch/csm6hg/genomes/sco_aug28/genes/${line}.fa

    mafft --auto \
    --thread 4 \
    /scratch/csm6hg/genomes/sco_aug28/genes/${line}.fa > \
    /scratch/csm6hg/genomes/sco_aug28/align/${line}.align.fa

done < /scratch/csm6hg/genomes/sco_aug28/sco.aug28.txt

# Align each SCO with mafft
while read -r line;
do
    # line=100070at6656
    echo ${line}

    singularity run /home/csm6hg/clipkit_latest.sif \
    clipkit -m smart-gap \
    /scratch/csm6hg/genomes/sco_aug28/align/${line}.align.fa

done < /scratch/csm6hg/genomes/sco_aug28/sco.aug28.txt

# Realignment of trimmed SCO with mafft
while read -r line;
do
    # line=100070at6656
    echo ${line}

    mafft --auto \
    --thread 4 \
    /scratch/csm6hg/genomes/sco_aug28/align/${line}.align.fa.clipkit > \
    /scratch/csm6hg/genomes/sco_aug28/align/${line}.realign.clip.fa 

done < /scratch/csm6hg/genomes/sco_aug28/sco.aug28.txt
