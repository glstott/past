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
    conda 'envs/simulation.yml'
    // Process for generating true trees given a number of taxa.

    input:
    val taxa
    val seed
    val out_tree
    val out_hist

    output:
    file "${out_tree}" into true_tree
    file "${out_hist}" into true_hist


    script:
    '''
    Rscript script.R -vanilla ${out_tree} ${out_hist} ${taxa} ${seed > 0 ? "${seed}" : ""}
    '''
}

process simulateSequences{
    conda "-c bioconda iqtree"
    // Process for simulating sequences given a true tree.

    input:
    file true_tree from true_tree
    val seed

    output:
    file "${true_tree.simpleName}.fasta" into sim_seqs
    //file "${out_meta}" into sim_meta

    script:
    '''
    iqtree2 --alisim ${true_tree.simpleName}.fasta -m JC --length 1500 -t ${true_tree} --out-format fasta  ${seed > 0 ? "-seed ${seed}" : ""}
    '''
}

workflow{
    // Define the workflow for the simulation.

    // Define the parameters for the simulation.
    taxa = 1000
    seed = 0
    out_tree = "true_tree.tre"
    out_hist = "true_hist.txt"
    //out_meta = "sim_meta.csv"

    // Run the processes for the simulation.
    trueTree(taxa, seed, out_tree, out_hist)
    simulateSequences(true_tree, out_meta)
}
```