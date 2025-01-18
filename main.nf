#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---------------------------
// (A) Define Channels
// ---------------------------

proteins = Channel
    .fromPath(params.protein_list)
    .splitText()
    .map { it.trim() }

// ---------------------------
// (B) Define Processes
// ---------------------------

// Process to extract longest transcript
process extract_longest_transcript {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/longest_orf", mode: 'copy'

    input:
        path input_faa

    output:
        path "primary_transcripts/*faa", emit: longest_transcripts

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate base

        # Extract longest transcript
        python ${params.OF_dir}/primary_transcript.py \\
            ${input_faa}
        """
}

// Define the workflow
workflow {
    // Extract longest transcript
    def orf = extract_longest_transcript(proteins)
}
