#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Run busco on proteomes 
process runBusco {

    input:
        path '*.faa'

    output:
        path "${sample_id}.busco"

    script:
        """
        cd /scratch/csm6hg/genomes/proteins_species/primary_transcripts/
        samp=${sample_id}
        species=$( echo \$samp | sed 's/.protein.faa//g' )
        echo \$species

        singularity run /home/csm6hg/busco_v5.4.7_cv1.sif \\
        busco \\
        -i /scratch/csm6hg/genomes/proteins_species/primary_transcripts/\${samp} \\
        -c 5 \\
        --out_path /scratch/csm6hg/genomes/proteins_species/primary_transcripts/ \\
        -l arthropoda_odb10 \\
        -o \${species} \\
        -m proteins
        """
}

// Plot BUSCO results
process plotBusco {
   
    input:
        path busco_results

    output:
        path "busco_plot.png"

    script:
        """
        module load singularity
        singularity run /home/csm6hg/busco_v5.4.7_cv1.sif \\
        python3 /scratch/csm6hg/scripts/generate_plot.py \\
        -wd /scratch/csm6hg/genomes/proteins_species/primary_transcripts/summary
        """
}

workflow {
    runBusco()
    plotBusco(runBusco.out.collect())
}
