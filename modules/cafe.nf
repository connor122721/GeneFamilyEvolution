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
    memory = '30 GB' // This is a memory intensive step
    time '1d'

    input:
        path(hog_input)
        path(tree)
        val(cafe_settings)

    output:
        path("r1*/*")

    script:
        """
        # Load modules
        module load apptainer

        cp ${hog_input} hog_gene_counts1.tsv 
        cp ${tree} mcmctree_busco_daphnia1.nwk 
        
        # Orthogroups from OrthoFinder
        apptainer run ${params.sif_dir}/cafe5_20240415.sif \\
            cafe5 \\
            -i hog_gene_counts1.tsv \\
            -t mcmctree_busco_daphnia1.nwk \\
            --cores 16 \\
            ${cafe_settings}
        
        # Extract significant gene families
        echo \$'#nexus\nbegin trees;' > Significant_trees.tre
        grep "*" */*asr.tre >> Significant_trees.tre
        echo "end;">> Significant_trees.tre
        awk '\$2 < .05 {print \$0}' */*family_results.txt > Sig_at_p.05.txt
        mv Sig* r1*/
        """
}

// GO term enrichment 
process enrichmentGO {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/go", mode: 'copy'

    input:
        path(cafe_out)

    output:
        path("*")

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1
        Rscript ${params.scripts_dir}/analyze_cafe_expand_contract_GO.R
        """
}
