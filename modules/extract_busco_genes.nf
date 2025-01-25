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
        path "*faa"
        path "*txt", emit: sco_list

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
        path sco_list

    output:
        path "genes/*"
        path sco_list, emit: sco

    script:
        """
        mkdir genes

        # Extract each SCO 
        while read -r line; do
            echo "\${line}"
            
            # Extract sequence
            ~/seqkit faidx ${params.out}/sco/${species}.rename.sco.faa "\${line}" > ./genes/"\${line}".${species}.faa

            # Rename
            sed -i "s/>"\${line}"/>${species}|"\${line}"/g" ./genes/"\${line}".${species}.faa

            # Add species header if file is empty (missing gene handler)
            # if [ ! -s ./genes/"\${line}".${species}.faa; then
            #     echo ">${species}|\${line}" > ./genes/"\${line}".${species}.faa
            # fi
        
        done < ${sco_list}
        """
}

// Aligning BUSCO genes
process alignBusco {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/sco", mode: 'copy'

    input:
        path(busco_file)

    output:
        path("align/*align.faa")

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 mafft/7.505

        mkdir align
        echo "Processing busco_file: ${busco_file}"

        # Align each SCO with mafft
        while read -r line; do
            echo "Aligning gene: \${line}"
            cat ${params.out}/sco/genes/"\${line}"*.faa > "\${line}".faa
            mafft --auto \\
               --thread ${params.threads} \\
               "\${line}".faa > \\
               align/"\${line}".align.faa
        done < ${busco_file}
        """
}

// Trims alignment of poor mapping
process trimAlign {
    
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/sco/trim", mode: 'copy'
    cpus = 1
    memory = '1GB'
    errorStrategy = 'finish'

    input:
        path(align_faa)

    output:
        tuple path("*.realign.clip.faa"),
            val("${params.out}/sco/sco.txt")

    script:
        """
        module load apptainer/1.3.4
        gene=\$( echo ${align_faa} | cut -f1 -d"." )

        # Align each SCO with mafft
        apptainer run ${params.sif_dir}/clipkit_latest.sif \\
            clipkit -m smart-gap \\
            ${align_faa} \\
            -o "\${gene}".clip.faa

        module load gcc/11.4.0 openmpi/4.1.4 mafft/7.505

        # Realignment of trimmed SCO with mafft
        mafft --auto \\
            "\${gene}".clip.faa > \\
            "\${gene}".realign.clip.faa
        """
}

// Make BUSCO trees
process runBuscoTrees {
    
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/trees", mode: 'copy'
    cpus = 1
    memory = '1GB'
    errorStrategy = 'finish'

    input:
        path(input_trim_faa)

    output:
        path("*.treefile")

    script:
        """
        gene=\$( echo ${input_trim_faa} | cut -f1 -d"." )

        # Run iqtree2
        ${params.iqtree_dir}/iqtree2 \\
            -redo \\
            -bb 1000 \\
            -s ${input_trim_faa} \\
            -T 1 \\
            --prefix \${gene}
        #cat trees/*.tre.treefile > trees/busco_sco_genes.tre
        """
}

// Add modules
include { makeConsensusMCMC } from './makeSpeciesTree.nf'

// Define the workflow
workflow {
        
    // Extract BUSCO genes
    def busco_genes_ext = extractBuscoGenes()

    Channel
        busco_genes_ext.sco_list
        .unique()
        .collect()
        .map { list -> list[0] }
        .set { sco_busco_results }

    // Extract and align BUSCO genes
    def ext_busco = extractBusco(species, 
                                 sco_busco_results)

    // Collect unique BUSCO result directories
    Channel
        ext_busco.sco
        .unique()
        .collect()
        .map { list -> list[0] }
        .set { ext_busco_results }

    //ext_busco_results.view()

    // Align BUSCO genes
    alignBusco(ext_busco_results)
    
    // Collect unique align's output
    //Channel align_busco.sco_align.unique().collect().map { list -> list[1] }
    Channel
        alignBusco.out
        .map { tuple -> tuple[0] }
        .flatten()    
        .set { align_busco_dir }

    align_busco_dir.view()

    // Trim and realign BUSCO genes
    trimAlign(align_busco_dir)

    // Each SCO
    Channel
       trimAlign.out
        .map { tuple -> tuple[0] }
        .flatten()
        .set { trim_faa }
    
    //trim_faa.view()

    // Make Busco trees
    runBuscoTrees(trim_faa)

    Channel
        runBuscoTrees.out
        .unique()
        .collect()
        .set { busco_trees }

    // Make consensus species-tree
    makeConsensusMCMC(busco_trees)
}