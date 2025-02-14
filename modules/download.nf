// Process to download sequences
process download_NCBI {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/genomes", mode: 'copy'
    errorStrategy 'ignore'

    input:
        path(ids_file)

    output:
        path("proteomes/*"), emit: protein_faa
        path("gtf/*")
        path("cds/*")

    script:
    """
    module load miniforge/24.3.0-py3.11
    source activate msprime_env

    python3 ${params.scripts_dir}/download_genomes.py \\
        --assembly_list ${ids_file} \\
        --out_dir downloads

    # Optional genome additions
    cp ${params.rawgenome}/pulex_euro/pulex.euro.protein.faa downloads/pulexeuro/pulexeuro.protein.faa
    cp ${params.rawgenome}/pulex_euro/Daphnia.aed.0.6.gff downloads/pulexeuro/pulexeuro.gff
    cp ${params.rawgenome}/pulex_euro/pulex.euro_cds_from_genomic_rename.fna downloads/pulexeuro/pulexeuro.cds_from_genomic.fna
    cp /project/berglandlab/daphnia_ref/Daphnia_annotation_PANTHER.xls downloads/pulexeuro/Daphnia_annotation_PANTHER.xls

    # Remove intermediate files
    rm downloads/*/*.gz

    # Get everything in the correct directories
    mkdir proteomes
    mkdir gtf
    mkdir cds

    mv downloads/*/*protein.faa proteomes/
    mv downloads/*/*gtf gtf/
    mv downloads/*/*gff gtf/
    mv downloads/*/*xls gtf/
    mv downloads/*/*cds_from_genomic.fna cds/
    mv downloads/*/*genomic.fna cds/
    """
}
