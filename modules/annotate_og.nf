#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Annotate Orthogroups
process annotateOrthogroups {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/orthogroups_annotation", mode: 'copy'

    input:
        path(fasta_dir)

    output:
        path("hogs.*.function.tsv")

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate msprime_env

        # Annotate Orthogroups
        annotate_orthogroups \\
            --orthogroups_tsv ${params.out}/longest_orf/primary_transcripts/OrthoFinder/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv \\
            --hog True \\
            --fasta_dir ${fasta_dir} \\
            --file_endings faa \\
            --out hogs.function.tsv \\
            --simple True
        """
}

 
// Annotate Orthogroups
process annotateGO {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/orthogroups_annotation", mode: 'copy'

    input:
        val(species)

    output:
        path("${species}.genes_orf_Pfam-A.tblout")
        path("${species}.genes_orf_Pfam-A.map")
        path("${species}.genes_orf_Pfam-A.map.gaf")
        path("*GOterm_mapping.tsv")

    script:
        """
        module load apptainer/1.3.4

        # Annotate genes with Pfam
        apptainer run ${params.sif_dir}/hmmer2go_latest.sif \\
            hmmer2go run \\
            -i ${params.out}/longest_orf/primary_transcripts/${species}.protein.faa \\
            -d ${params.pfamdb_dir}/Pfam-A.hmm \\
            -n ${params.threads} \\
            -o ${species}.genes_orf_Pfam-A.tblout

        ## Find Pfam -> GO term mapping
        apptainer run ${params.sif_dir}/hmmer2go_latest.sif \\
            hmmer2go mapterms \\
            -i ${species}.genes_orf_Pfam-A.tblout \\
            -o ${species}.genes_orf_Pfam-A.map \\
            --map

        ## Create GAF file
        apptainer run ${params.sif_dir}/hmmer2go_latest.sif \\
            hmmer2go map2gaf \\
            -i ${species}.genes_orf_Pfam-A.map \\
            -o ${species}.genes_orf_Pfam-A.map.gaf \\
            -s ${species}
        """
}



