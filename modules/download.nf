// Process to download sequences
process download_sequences {
    input:
    path ids_file

    output:
    path "sequences/*.faa" into downloaded_sequences
    path "annotations/*.gtf" into downloaded_gtf
    path "annotations/*.gff" into downloaded_gff

    script:
    """
    mkdir -p sequences annotations
    while read id; do
        # Fetch protein sequences
        efetch -db protein -id \$id -format fasta > sequences/\$id.faa
        # Fetch GTF and GFF annotations
        efetch -db nuccore -id \$id -format gtf > annotations/\$id.gtf
        efetch -db nuccore -id \$id -format gff > annotations/\$id.gff
    done < \$ids_file
    """
}

// Process to extract longest transcript
process extract_longest_transcript {
    input:
    path downloaded_sequences

    output:
    path "longest_transcripts/*.faa"

    script:
    """
    mkdir -p longest_transcripts
    for f in sequences/*.faa; do
        python ~/OrthoFinder_source/tools/primary_transcript.py \$f
    done
    """
}

// Define the workflow
workflow {
    // Download sequences
    download_sequences(params.ids_file)

    // Extract longest transcript
    extract_longest_transcript()
}
