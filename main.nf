#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---------------------------
// (A) Define Channels
// ---------------------------

// Extract species names from protein file paths
Channel
    .fromPath(params.protein_list)
    .splitText()
    .map { file -> 
        def species = file.tokenize('/').last().tokenize('.').first()
        tuple(file, species)
    }
    .set { proteomes }

// ---------------------------
// (B) Define Processes
// ---------------------------

// Process to extract longest transcript
process extract_longest_transcript {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/longest_orf", mode: 'copy'

    input:
        tuple path(input_faa), 
        val(species)

    output:
        val "${params.out}/longest_orf/primary_transcripts/", emit: output_dir
        path "primary_transcripts/${species}*faa", emit: output_faa

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate base

        # Extract longest transcript
        python ${params.OF_dir}/primary_transcript.py \\
            ${input_faa}
        """
}

// Module imports
include { runBusco; plotBusco } from './modules/run_busco.nf'
include { runOrthoFinder } from './modules/run_orthofinder.nf'

// Define the workflow
workflow {
        
    // Extract longest transcript
    def orf = extract_longest_transcript(proteomes)

    // Run BUSCO on proteomes
    def busco = runBusco(proteomes)

    // Collect unique BUSCO result directories
    busco.busco_dir
        .collect()
        .unique()
        .set { busco_results }
    
    // Plot BUSCO results after all BUSCO processes are complete
    plotBusco(busco_results)

    // Collect all longest transcripts paths into a list after all extract_longest_transcript processes are complete
    orf.output_dir
        .unique()
        .set { longest_transcripts_dirs }

    // Run OrthoFinder on longest transcripts
    runOrthoFinder(longest_transcripts_dirs)
}
