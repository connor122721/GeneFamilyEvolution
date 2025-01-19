#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Extracting BUSCO genes
process extractBuscoGenes {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/sco", mode: 'copy'

    input:
        path busco_files

    output:
        path "busco_output"

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        Rscript ${params.scripts_dir}/extract_busco_genes.R
        """
}
