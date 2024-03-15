#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Authors: Guppy Stott
// Purpose: Simulate a true tree and sequences for testing subsampling methods

// Print some log information for the user.
log.info """
Discrete Phylogeography Simulation
==================================
Project : $workflow.projectDir
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
==================================
"""

process trueTree{
    //conda "-c bioconda -c conda-forge r-base r-ape=5.7 r-diversitree r-phangorn=2.11 bioconductor-treeio=1.26.0"
    conda "$workflow.projectDir/envs/simulation.yml"
    // Process for generating true trees given a number of taxa.

    input:
    val taxa
    val seed
    val prefix

    output:
    path "${prefix}-${seed}-${taxa}.nwk", emit: ttree


    script:
    """
    pwd
    Rscript $workflow.projectDir/scripts/simulate.r --vanilla ${prefix}-${seed}-${taxa}.nwk ${taxa} ${seed}
    """
}

process simulateSequences{
    conda "-c bioconda iqtree"
    // Process for simulating sequences given a true tree.

    input:
    path true_tree
    val seed

    output:
    path "${true_tree.simpleName}.fa"
    //file "${out_meta}" into sim_meta

    script:
    """
    iqtree2 --alisim ${true_tree.simpleName} -m HKY --length 1500 -t ${true_tree} --out-format fasta  ${seed}
    """
}

process simulateConvenienceSampling {
    conda "-c bioconda seqkit"
    //simulate the acquisition of a convenience sample from the simulated sequences
    input:
    path seqs
    val seed
    val prefix
    val p

    output:
    path "${prefix}-${seed}-convenience.fa"

    script:
    """
    seqkit sample -p ${p} -s ${seed} ${seqs} > ${prefix}-${seed}-convenience.fa
    """
}

process simulateBiasedSampling {
    conda "-c bioconda seqkit"
    // separate fasta into files by the second delimiter, then sample from each file
    input:
    path seqs
    val seed
    val prefix
    val p

    output:
    path "${prefix}-${seed}-biased.fa"

    script:
    """
    seqkit split ${seqs} -i --id-regexp "^[^_]*_([^_]*)_"
    for file in *part_.fa; do
        seqkit sample -p ${p} -s ${seed} $file >> ${prefix}-${seed}-biased.fa
    done
    """
}

workflow{
    // Define the workflow for the simulation.

    // Define the parameters for the simulation.
    taxa = 100
    seed = 12
    prefix = "test0315"
    //out_meta = "sim_meta.csv"

    // Run the processes for the simulation.
    trueTree(taxa, seed, prefix)
    simulateSequences(trueTree.out.ttree, seed)
    simulateConvenienceSampling(simulateSequences.out, seed, prefix, 0.1)
    simulateBiasedSampling(simulateSequences.out, seed, prefix, 0.1)
}
