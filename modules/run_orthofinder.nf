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

        # Run OrthoFinder
        orthofinder \\
            -f ${longest_transcripts_dir} \\
            -t ${params.of_threads}
        """
}
