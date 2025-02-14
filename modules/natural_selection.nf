#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Translate and detect selection 
process runSelection {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/selection", mode: 'copy'
    errorStrategy 'finish'

    input:
        path(gene_fa)

    output:
        path("*result*")

    script:
        """
        # Extract gene id from the input file name
        id=\$(basename ${gene_fa} .fa)
        echo "Processing gene file: ${gene_fa} ; Orthogroup = \${id}"

        #id="OG0003402"

        # Genomes
        genome="${params.rawgenome}"

        # Create and move to working directory
        if [ ! -d mega_OG/\${id} ]; then
            mkdir -p mega_OG/\${id}
        fi
        cd mega_OG/\${id}

        ################ Arrange codon MSAs ###################
        # Make list of genes from the input gene FASTA
        cat ../../${gene_fa} | grep ">" | sed 's/>//g' > \${id}.genes

        # Extract gene sequences
        for i in \${genome}/*/*cds_from_genomic_rename.fna; do
            echo \${i}
            # Extract gene sequences (DNA)
            ~/seqkit faidx \${i} -l \${id}.genes >> \${id}.fa
        done

        # Align the nuc
        module load gcc/11.4.0 openmpi/4.1.4 mafft/7.505
        mafft --auto \\
            \${id}.fa > \\
            \${id}.align.fa

        # Align aa
        mafft --auto \\
            ../../${gene_fa} > \\
            \${id}.align.faa

        # Download pal2nal: http://www.bork.embl.de/pal2nal/#RunP2N        
        module load miniforge/24.3.0-py3.11
        source activate papaml

        # Translate and make input for codeml
        ~/phd_software/pal2nal.pl \\
            \${id}.align.faa \\
            \${id}.align.fa \\
            -output fasta -nogap -nomismatch > \\
            \${id}.trans.fa

        # Clip out poorly aligned regions
        module load apptainer
        apptainer run ${params.sif_dir}/gblocks_latest.sif Gblocks \\
            \${id}.trans.fa \\
            -t=c \\
            -p=n \\
            -b5=h \\
            -e=.gb

        # Make fasta > 1 line per sequence
        perl ${params.scripts_dir}/one_line_fasta.pl \${id}.trans.fa.gb

        # Remove spaces introduced through Gblocks processing
        cat \${id}_one_line.fa | tr -d " " > \${id}.trans.clip.fa

        # Run iqtree2 for a gene tree
        ${params.iqtree_dir}/iqtree2 \\
            -redo \\
            -s \${id}.trans.clip.fa \\
            -T 1 \\
            --prefix \${id}.tre

        # Remove dup sequences
        hyphy ~/phd_software/hyphy-analyses/remove-duplicates/remove-duplicates.bf \\
            --msa \${id}.trans.clip.fa \\
            --tree \${id}.tre.treefile \\
            --output \${id}.uniques.nxh

        # Conditional run Hyphy
        # Run aBSREL: install hyphy -> conda install -c bioconda hyphy
        if [ -f \${id}.uniques.nxh ]; then
            echo "Duplicates found, running hyphy."
            hyphy absrel \\
                --alignment \${id}.uniques.nxh \\
                --tree \${id}.tre.treefile \\
                --output \${id}.ABSREL.json
        else
            echo "No duplicates found."
            hyphy absrel \\
                --alignment \${id}.trans.clip.fa \\
                --tree \${id}.tre.treefile \\
                --output \${id}.ABSREL.json
        fi
        """
}

// Workflow: pipe every gene FASTA file from the orthofinder output to a runSelection job
workflow {

    // Create a channel from gene FASTA files 
    Channel.fromPath("${params.out}/longest_orf/primary_transcripts/OrthoFinder/Results_Feb04/Orthogroup_Sequences/*.fa")
           .set { geneChannel }

    runSelection(geneChannel)
    
}