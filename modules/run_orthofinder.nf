#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Run OrthoFinder on longest transcripts
process runOrthoFinder {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/orthofinder", mode: 'copy'

    input:
        path longest_transcripts_dir

    output:
        path "*"

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate msprime_env

        # Exception species conditional statement
        if [[ ${longest_transcripts_dir}/pulexeuro_orf.protein.faa ]]; then
            mv ${longest_transcripts_dir}/pulexeuro_orf.protein.faa \\
            ${longest_transcripts_dir}/pulexeuro.protein.faa
        fi

        # Run OrthoFinder
        orthofinder \\
            -f ${longest_transcripts_dir} \\
            -t ${params.of_threads}
        """
}
