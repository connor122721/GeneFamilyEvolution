#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Filter HOGs
process filterHogs {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/cafe", mode: 'copy'

    input:
        path(hog_input)

    output:
        path("hog_gene_counts.tsv")

    script:
        """
        # Run prep script
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1
        Rscript ${params.scripts_dir}/FilterOrthofinderHogs.R
        """
}

// Run OrthoFinder on longest transcripts
process runCAFE {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/cafe", mode: 'copy'
    cpus 16
    memory = '50 GB' // This is a memory intensive step
    time '1d' // This also takes awhile to finish

    input:
        path(hog_input)
        path(spp_tree)

    output:
        path("*")

    script:
        """
        # Load modules
        module load gcc intel

        # Cafe5 executable
        cafe=".local/bin/cafe5"

        # Ultrametric time-calibrated tree from MCMCtree
        tree="/scratch/csm6hg/cafe_aug28/mcmctree_busco_daphnia.nwk"

        # Orthogroups from OrthoFinder
        ${cafe} \\
            -i ${ortho} \\
            -t ${tree} \\
            -o r1_hogs \\
            --cores 16

        ${cafe} \\
            -i ${ortho} \\
            -t ${tree} \\
            -o r1_k3_hogs \\
            -k 3 \\
            --cores 16

        ${cafe} \\
            -i ${ortho} \\
            -t ${tree} \\
            -o r1_k3_p_hogs \\
            -k 3 \\
            -p \\
            --cores 16 

        ${cafe} \\
            -i ${ortho} \\
            -t ${tree} \\
            -o r1_e_hogs \\
            -e \\
            --cores 16
        """
}
