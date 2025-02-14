// GO term enrichment 
process enrichmentGO {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/gene_family_GO", mode: 'copy'
    cpus = 4
    memory = '15GB'

    input:
        path(cafe_out)

    output:
        path("*.rds")
        path("*.txt")
        path("*.pdf")

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1
        Rscript ${params.scripts_dir}/analyze_cafe_expand_contract_GO.R \\
        --cafe_dir ${cafe_out}
        """
}

workflow {

    // Create a channel from CAFE
    Channel.fromPath("${params.out}/cafe")
           .set { cafe_ch }

    enrichmentGO(cafe_ch)
    
}