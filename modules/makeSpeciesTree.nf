// Make consensus species tree and prep for MCMCtree
process makeConsensusMCMC {
    
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/spp_tree", mode: 'copy'

    input:
        path(sco_tree_list)

    output:
        path("*.nwk")
        path("busco_sco_genes.tre")
        path("sco_species_astralpro")

    script:
        """
        # Make full gene tree
        cat *treefile > busco_sco_genes.tre

        # Get rid of sco geneid
        while read -r line; do
            echo \${line}
            sed -i 's/'\${line}'//g' busco_sco_genes.tre
        done < *treefile

        # Remove "|" character
        sed -i 's/|:/:/g' busco_sco_genes.tre
        sed -i 's/|2:/:/g' busco_sco_genes.tre

        # Run Astral-pro to create species tree
        ~/phd_software/ASTER-1.15/bin/astral-pro \\
            -i busco_sco_genes.tre \\
            -R \\
            --root "melanogaster" \\
            -o sco_species_astralpro

        # Make prep tree for MCMCtree
        module load miniforge/24.3.0-py3.11
        source activate msprime_env

        # Add time constraints
        python ${params.scripts_dir}/mcmctree_prep.py \\
            --left_species pulexnam \\
            --right_species magna \\
            --lower_bound 130 \\
            --upper_bound 150 \\
            --tree sco_species_astralpro \\
            | python ${params.scripts_dir}/mcmctree_prep.py \\
            --left_species magna \\
            --right_species sinensis \\
            --lower_bound 21.5 \\
            --upper_bound 22.4 \\
            --tree - \\
            | python ${params.scripts_dir}/mcmctree_prep.py \\
            --left_species magna \\
            --right_species melanogaster \\
            --lower_bound 474.8 \\
            --upper_bound 530 \\
            --tree - \\
            | python ${params.scripts_dir}/mcmctree_prep.py \\
            --left_species anopheles \\
            --right_species melanogaster \\
            --lower_bound 216.6 \\
            --upper_bound 265.7 \\
            --tree - \\
            --add_header \\
            > sco_daphnia_mcmc_1.25.25.nwk
        """
}