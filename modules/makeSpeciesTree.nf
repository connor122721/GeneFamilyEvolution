// Make consensus species tree and prep for MCMCtree
process makeConsensusMCMC {
    
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/spp_tree", mode: 'copy'

    input:
        path(sco_tree_list)

    output:
        path("*.nwk"), emit: nwk
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
        done < ${params.sco_data}

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
        source activate base

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
            --add_header \\
            > sco_daphnia_mcmc_1.25.25.nwk
        """
}

// Prep MCMCtree input
process prepMCMCtree {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/spp_tree", mode: 'copy'

    input:
        path(nwk)

    output:
        path("sco_daphnia_mega_clip.faa"), emit: ali

    script:
        """
        # Make mega sequence alignment for MCMCtree
        while read -r gene; do
            echo \${gene}
                    
            # Rename
            sed "s/|\${gene}//g" \\
                ${params.out}/sco/trim/\${gene}.realign.clip.faa > \\
                \${gene}.align.clip.faa

        done < ${params.out}/sco/sco.txt

        # Concatenate aligned genes
        ~/seqkit concat *.align.clip.faa \\
            -o sco_daphnia_mega_clip.faa
        """
}

// Run MCMCtree
process MCMCTREE {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/mcmctree", mode: 'copy'
    cpus 8
    memory = '50 GB' // This is a memory intensive step
    time '1d' // This also takes awhile to finish

    input:
        path(ali)
        path(nwk)
        val(replicate)

    output:
        path("FigTree_*")
        path("mcmctreeout*")
        path("mcmc_*.txt")

    script:
        """
        # Run PAML package MCMCtree
        export PATH=~/phd_software/paml-4.10.7/bin/:$PATH

        cp ${params.wd}/bin/mcmctree.ctl .
        sed -i "s/outfile = mcmctreeout/outfile = mcmctreeout_${replicate}/g" mcmctree.ctl

        # Run mcmctree - first run to estimate in.BV
        mcmctree mcmctree.ctl

        # modify control file
        mv out.BV in.BV
        sed -i "s/usedata = 3/usedata = 2/g" mcmctree.ctl

        # Rerun to get divergence times
        mcmctree mcmctree.ctl
        
        # Rename for replicates
        mv FigTree.tre FigTree_${replicate}.tre
        mv mcmc.txt mcmc_${replicate}.txt
        """
}

// Plot MCMCtree output
process plotMCMCtree {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/mcmctree", mode: 'copy'
    memory = '30 GB' // This is a memory intensive step

    input:
        path(Figtree)
        path(mcmctree_log)
        path(mcmctree_output)

    output:
        path("MCMCtree.pdf")
        path("mcmc.full.pdf")
        path("mcmctree_busco_daphnia.nwk"), emit: cafe_input_tree

    script:
        """
        # Run analysis script
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1
        Rscript ${params.scripts_dir}/processMCMCtree.R
        mv Rplots.pdf MCMCtree.pdf
        
        """
}
