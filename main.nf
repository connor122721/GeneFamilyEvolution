#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---------------------------
// (A) Define Processes
// ---------------------------

// Process to extract longest transcript
process extract_longest_transcript {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/longest_orf", mode: 'copy'

    input:
        tuple path(input_faa), 
            val(species)

    output:
        val("${params.out}/longest_orf/primary_transcripts/"), emit: output_dir
        path("primary_transcripts/*protein.faa"), emit: output_faa
        val("${species}"), emit: species

    script:
        """
        # Extracting the longest ORF for expection species (European Daphnia pulex)
        # This is due to the reference genome not being annotated
        if [ ${species} == "pulexeuro" ]; then
            module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1
            module load samtools
            mkdir primary_transcripts
            cd primary_transcripts
            samtools faidx ../${input_faa}
            cat ../${input_faa}.fai | cut -f1 > pulexeuro.sequences
            Rscript ${params.scripts_dir}/longest_orf_europulex.R
            cp "${input_faa}" .
            rm pulexeuro_orf.protein.faa
        else
            # Extract longest transcript
            module load miniforge/24.3.0-py3.11
            source activate base
            python ${params.OF_dir}/primary_transcript.py \\
                ${input_faa}
        fi 
        """
}

// Module imports
include { download_NCBI } from './modules/download.nf'
include { runBusco; plotBusco } from './modules/run_busco.nf'
include { runOrthoFinder } from './modules/run_orthofinder.nf'
include { annotateGO } from './modules/annotate_og.nf'

// Define the workflow
workflow {

    // Run download workflow which uses params.ncbi_genomes
    def downloaded_files = download_NCBI(params.species_list)

    // Map downloaded file paths to tuple (file, species)
    def proteomes = downloaded_files.protein_faa
        .flatten()
        .map { file ->
            def filename = file.toString().trim().tokenize('/').last()
            def species = filename.tokenize('.')[0]  // species name is the first token before the period
            tuple(file, species)
        }

    // Extract longest transcript
    def orf = extract_longest_transcript(proteomes)
   
    // Run BUSCO on proteomes
    def busco = runBusco(proteomes)

     // Annotate GO terms
    def annot_go = annotateGO(busco.species)

    // Collect unique BUSCO result directories
    def busco_results = busco.busco_dir.unique()
        .collect()
        .map { list -> list[0] }
    
    // Plot BUSCO results after all BUSCO processes are complete
    plotBusco(busco_results)

    // Collect all longest transcripts paths into a list after all extract_longest_transcript processes are complete
    def longest_transcripts_dirs = orf.output_dir
        .unique()
        .collect()
        .map { list -> list[0] }

    // Run OrthoFinder on longest transcripts
    def orthofinder = runOrthoFinder(longest_transcripts_dirs)
}
