#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Translate and detect selection 
process runSelection {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/selection", mode: 'copy'

    input:
        val(of_dir)

    output:
        path("*result*")

    script:
        """
        # Conditional make mega list for codeml
        if [ -f mega.list.OG ]; then
            echo "Mega list is there."
        else
            echo "Empty mega list."
            ls -d ${of_dir}/*fa > \\
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

        ################ Arrange codon MSAs ###################

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
        mafft --auto \\
            ${id}.fa > \\
            ${id}.align.fa

        # Align aa
        mafft --auto \\
            ${gene} > \\
            ${id}.align.faa

        # Download pal2nal: http://www.bork.embl.de/pal2nal/#RunP2N
        # Translate and make input for codeml
        ~/phd_software/pal2nal.pl \\
            ${id}.align.faa \\
            ${id}.align.fa \\
            -output fasta -nogap -nomismatch > \\
            ${id}.trans.fa

        # Clip out poorly aligned regions
        module load miniforge/24.3.0-py3.11
        source activate papaml

        module load apptainer
        apptainer run ~/gblocks_latest.sif Gblocks \\
            ${id}.trans.fa \\
            -t=c \\
            -p=n \\
            -b5=h \\
            -e=.gb

        # Make fasta > 1 line per sequence
        perl ${params.scripts_dir}/one_line_fasta.pl ${id}.trans.fa.gb

        # Remove spaces introduced through Gblocks processing
        cat ${id}_one_line.fa | tr -d " " > ${id}.trans.clip.fa

        # Run iqtree2 for a gene tree
        ${params.iqtree_dir}/iqtree2 \\
            -redo \\
            -s ${id}.trans.clip.fa \\
            -T 1 \\
            --prefix ${id}.tre

        # Remove dup sequences
        hyphy ~/phd_software/hyphy-analyses/remove-duplicates/remove-duplicates.bf \\
            --msa ${id}.trans.clip.fa \\
            --tree "/scratch/csm6hg/codeml/mega_OG/${id}/${id}.tre.treefile" \\
            --output ${id}.uniques.nxh

        # Conditional run Hyphy
        # Run aBSREL: install hyphy -> conda install -c bioconda hyphy
        if [ -f ${id}.uniques.nxh ]; then
            echo "Duplicates found, running hyphy."
            hyphy absrel \\
                --alignment ${id}.uniques.nxh \\
                --tree ${id}.tre.treefile \\
                --output ${id}.ABSREL.json
        else
            echo "No duplicates found."
            hyphy absrel \\
                --alignment ${id}.trans.clip.fa \\
                --tree ${id}.tre.treefile \\
                --output ${id}.ABSREL.json
        fi
        """
}