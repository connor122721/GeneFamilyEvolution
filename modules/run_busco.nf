#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Download Busco dataset
process downloadBusco {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/busco_dataset", mode: 'copy'

    input:
        val dataset

    output:
        path "*"

    script:
        """
        module load apptainer/1.3.4

        # Run BUSCO
        apptainer run ${params.sif_dir}/busco_v5.4.7_cv1.sif \\
            busco \\
            --download ${dataset}
        """
}

// Run busco on proteomes 
process runBusco {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/busco_results", mode: 'copy'

    input:
        tuple path(samp), 
        val(species)

    output:
        val "${params.out}/busco_results/", emit: busco_dir
        path "short_summary*txt"
        path "full_table_${species}.tsv"
        val "${species}", emit: species

    script:
        """
        module load apptainer/1.3.4

        # Determine BUSCO lineage
        #if [ ${species} == "caenorhabditis" ]; then
        #    lineage="nematoda_odb10"
        #else
        #    lineage="${params.busco_lineage}"
        #fi

        # Run BUSCO with local data
        apptainer run ${params.sif_dir}/busco_v5.4.7_cv1.sif \\
            busco \\
            -i ${params.out}/genomes/proteomes/${samp} \\
            -c ${params.threads} \\
            -l ${params.busco_lineage} \\
            -o ${species} \\
            -m proteins \\
            --offline \\
            --download_path ${params.busco_data}

        # Rename full_table.tsv to include species name
        mv ${species}/run_*/full_table.tsv full_table_${species}.tsv
        mv ${species}/short_summary.specific*${species}.txt .
        """
}

// Plot BUSCO results
process plotBusco {
   
    // Set the error strategy to ignore failures in this process, dumb plotting errors 
    errorStrategy 'ignore'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/busco_plot", mode: 'copy'

    input:
        path busco_results

    output:
        path "busco_plot.png"

    script:
        """
        module load apptainer/1.3.4

        # Run plot creation script
        apptainer run ${params.sif_dir}/busco_v5.4.7_cv1.sif \\
            python3 ${params.scripts_dir}/generate_plot.py \\
            -wd ${params.out}/${busco_results} \\
            --no_r
        
        # Run Busco script
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1
        Rscript ${params.out}/${busco_results}/busco_figure.R
        """
}
