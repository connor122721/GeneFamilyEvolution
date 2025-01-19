#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Run busco on proteomes 
process runBusco {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/busco_results", mode: 'copy'

    input:
        tuple path(samp), 
        val(species)

    output:
        val "${params.out}/busco_results/", emit: busco_dir
        path "*/short_summary*txt"
        path "*/run_arthropoda_odb10/full_table.tsv"

    script:
        """
        module load apptainer/1.3.4

        # Run BUSCO with local data
        apptainer run ${params.sif_dir}/busco_v5.4.7_cv1.sif \\
            busco \\
            -i ${params.rawpro}/${samp} \\
            -c ${params.threads} \\
            -l ${params.busco_lineage} \\
            -o ${species} \\
            -m proteins \\
            --offline \\
            --download_path ${params.busco_data}
        """
}

// Plot BUSCO results
process plotBusco {
   
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/busco_plot", mode: 'copy'

    input:
        path busco_results

    output:
        path "busco_plot.png"

    script:
        """
        module load apptainer/1.3.4

        apptainer run ${params.sif_dir}/busco_v5.4.7_cv1.sif \\
            python3 ${params.scripts_dir}/generate_plot.py \\
            -wd ${params.out}/${busco_results}
        """
}
