#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Run OrthoFinder on longest transcripts
process runOrthoFinder {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/orthofinder", mode: 'copy'
    cpus 16
    memory = '50 GB' // This is a memory intensive step
    time '1d' // This also takes awhile to finish

    input:
        path longest_transcripts_dir

    output:
        path "orthofinder_results*"
        path "hogs.function.tsv"
        val "${params.out}/longest_orf/primary_transcripts", emit: fasta_dir

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate msprime_env

        # Exception species conditional statement
        if [[ ${longest_transcripts_dir}/pulexeuro_rename.protein.faa ]]; then
            mv ${longest_transcripts_dir}/pulexeuro_rename.protein.faa \\
            ${longest_transcripts_dir}/pulexeuro.protein.faa
        fi

        # Run OrthoFinder
        orthofinder \\
            -f ${longest_transcripts_dir} \\
            -t ${params.of_threads} \\
            -a ${params.of_threads} \\
            -o orf

        mv orf/Results*/ \\
            orthofinder_results

        # Annotate Orthogroups
        annotate_orthogroups \\
            --orthogroups_tsv orthofinder_results/Phylogenetic_Hierarchical_Orthogroups/N0.tsv \\
            --hog True \\
            --fasta_dir ${longest_transcripts_dir} \\
            --file_endings faa \\
            --out hogs.function.tsv \\
            --simple True
        """
}
