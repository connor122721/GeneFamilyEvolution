#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---------------------------
// (A) Define Channels
// ---------------------------

// Extract species names from protein file paths
// Define the channel from the input file
Channel
    .fromPath(params.protein_list)
    .splitText()
    .map { line -> 
        def (file, group) = line.trim().tokenize('\t')
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
        path "primary_transcripts/${species}*.protein.faa", emit: output_faa
        val "${species}", emit: species

    script:
        """
        # Extracting the longest ORF for expection species 
        if [ ${species} == "pulexeuro" ]; then
            module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1
            mkdir primary_transcripts
            cd primary_transcripts
            Rscript ${params.scripts_dir}/longest_orf_europulex.R
            cp "/project/berglandlab/connor/genomes/proteins_species/${input_faa}" .
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
include { runBusco; plotBusco } from './modules/run_busco.nf'
include { runOrthoFinder } from './modules/run_orthofinder.nf'
include { annotateOrthogroups; annotateGO } from './modules/annotate_og.nf'
include { extractBuscoGenes } from './extract_busco_genes.nf'

// Define the workflow
workflow {
        
    // Extract longest transcript
    def orf = extract_longest_transcript(proteomes)

    // Annotate GO terms
    def annot_go = annotateGO(orf.species)

    // Run BUSCO on proteomes
    def busco = runBusco(proteomes)

    // Collect unique BUSCO result directories
    Channel
        busco.busco_dir
        .unique()
        .collect()
        .map { list -> list[0] }
        .set { busco_results }
    
    // Plot BUSCO results after all BUSCO processes are complete
    plotBusco(busco_results)

    // Collect all longest transcripts paths into a list after all extract_longest_transcript processes are complete
    Channel
        orf.output_dir
        .unique()
        .collect()
        .map { list -> list[0] }
        .set { longest_transcripts_dirs }

    // Run OrthoFinder on longest transcripts
    def orthofinder = runOrthoFinder(longest_transcripts_dirs)

    // Annotate Orthogroups and generate plots
    def annot_of = annotateOrthogroups(orthofinder)

}
