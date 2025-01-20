#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Extract species names from protein file paths
Channel
    .fromPath(params.protein_list)
    .splitText()
    .map { file -> 
        def species = file.tokenize('/').last().tokenize('.').first()
        species
    }
    .set { species }

// Extracting BUSCO genes
process extractBuscoGenes {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/sco", mode: 'copy'

    output:
        path "*sco.faa", emit: busco_faa
        path "*"

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        Rscript ${params.scripts_dir}/extract_busco_genes.R \\
            --max_missing_ingroup ${params.max_missing_ingroup} \\
            --max_missing_outgroup ${params.max_missing_outgroup} \\
            --metadata ${params.protein_list}
        """
}

// Extracting and aligning BUSCO genes
process extractBusco {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/sco", mode: 'copy'

    input:
        val species
        path busco_file

    output:
        path "genes/*"
        val "${params.out}/sco/sco.txt", emit: sco

    script:
        """
        mkdir genes

        sco_data=${params.sco_data}

        # Extract each SCO 
        while read -r line; do
            echo "\${line}"
            
            # Extract sequence
            ~/seqkit faidx ${params.out}/sco/${species}.rename.sco.faa "\${line}" > ./genes/"\${line}".${species}.faa

            # Rename
            sed -i "s/>"\${line}"/>${species}|"\${line}"/g" ./genes/"\${line}".${species}.faa
        done < "\${sco_data}"
        """
}

// Extracting and aligning BUSCO genes
process alignBusco {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/sco", mode: 'copy'

    input:
        path busco_file

    output:
        path "align/*"

    script:
        """
        module load mafft/7.505

        mkdir align

        # Align each SCO with mafft
        while read -r line; do
            echo "\${line}"

            cat ${params.out}/genes/"\${line}"*.faa > "\${line}".faa

            mafft --auto \\
            --thread ${params.threads} \\
            "\${line}".faa > \\
            ./align/"\${line}".align.faa

        done < ${params.sco_data}
        """
}

// Define the workflow
workflow {
        
    // Extract BUSCO genes
    def busco_genes = extractBuscoGenes()

    // Extract and align BUSCO genes
    def ext_busco = extractBusco(species, params.sco_data)

    // Collect unique BUSCO result directories
    ext_busco.sco
        .unique()
        .collect()
        .set { ext_busco_results }

    // Align BUSCO genes
    def align_busco = alignBusco(ext_busco_results)

}