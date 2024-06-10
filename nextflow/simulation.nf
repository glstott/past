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

process trueTree {
    //conda "-c bioconda -c conda-forge r-base r-ape=5.7 r-diversitree r-phangorn=2.11 bioconductor-treeio=1.26.0"
    conda "$workflow.projectDir/envs/simulation.yml"
    publishDir "$workflow.projectDir/../raw/", mode: 'copy'
    // Process for generating true trees given a number of taxa.

    input:
    val taxa
    val seed
    val prefix

    output:
    path "${prefix}-${seed}-${taxa}.nwk", emit: ttree
    path "${prefix}-${seed}-${taxa}.csv", emit: tmeta


    script:
    """
    pwd
    Rscript $workflow.projectDir/scripts/simulate.r --vanilla ${prefix}-${seed}-${taxa}.nwk ${prefix}-${seed}-${taxa}.csv ${taxa} ${seed}
    """
}

process simulateSequences {
    conda "-c bioconda iqtree"
    publishDir "$workflow.projectDir/../raw/", mode: 'copy'
    // Process for simulating sequences given a true tree.

    input:
    path true_tree
    val seed

    output:
    path "${true_tree.simpleName}.fasta", emit: simseq
    //file "${out_meta}" into sim_meta

    script:
    """
    iqtree2 --alisim ${true_tree.simpleName} -m HKY --length 1500 -t ${true_tree} --out-format fasta  ${seed}
    """
}

process splitSamples {
    conda "-c bioconda seqkit"
    // Process for splitting a fasta file into multiple files.

    input:
    path seqs

    output:
    path "*part_.fa", emit: parts

    script:
    """
    seqkit split ${seqs} -i --id-regexp "^[^_]*_([^_]*)_"
    """
}

process simulateSimpleSampling {
    conda "-c bioconda seqkit"
    publishDir "$workflow.projectDir/../raw/", mode: 'copy'
    //simulate the acquisition of a convenience sample from the simulated sequences
    input:
    path seqs
    val seed
    val prefix
    val p
    path metadata

    output:
    path "${prefix}-${seed}-convenience.fasta", emit: simple

    script:
    """
    seqkit sample -p ${p} -s ${seed} ${seqs} > ${prefix}-${seed}-convenience.fasta
    grep "^>" ${prefix}-${seed}-convenience.fasta | sed 's/>//' > headers.txt
    """
}

process simulateBiasedSampling {
    conda "-c bioconda seqkit"
    publishDir "$workflow.projectDir/../raw/", mode: 'copy'
    // separate fasta into files by the second delimiter, then sample from each file
    input:
    path seqs
    val seed
    val prefix
    val p

    output:
    path "*part_.fa", emit: parts
    path "${prefix}-${seed}-biased.fa", emit: biased

    script:
    """
    seqkit split ${seqs} -i --id-regexp "^[^_]*_([^_]*)_"
    i=0

    for file in *part_.fa; do
        seqkit sample -p ${p[i]} -s ${seed} $file >> ${prefix}-${seed}-biased.fa
        i += 1
    done
    """
}

process runSimpleSampling {
    publishDir "$workflow.projectDir/../sampled/", mode: 'copy'
    input:
    path seqs
    val seed
    val n
    val prefix
    path metadata

    output:
    path "${prefix}-${seed}-simple.fasta", emit: seq
    path "${prefix}-${seed}-simple.csv", emit: meta

    script:
    """
    Rscript $workflow.projectDir/scripts/simpleSample.R ${metadata} id Collection_Date ${seqs} ${prefix}-${seed}-simple.csv ${prefix}-${seed}-simple.fasta ${n} ${seed}
    """
}

process runStratifiedSampling {
    publishDir "$workflow.projectDir/../sampled/", mode: 'copy'
    input:
    path seqs
    val seed
    val n
    val prefix
    path metadata

    output:
    path "${prefix}-${seed}-stratified.fasta", emit: seq
    path "${prefix}-${seed}-stratified.csv", emit: meta

    script:
    """
    Rscript $workflow.projectDir/scripts/stratifiedSample.R ${metadata} id Collection_Date ${seqs} ${prefix}-${seed}-stratified.csv ${prefix}-${seed}-stratified.fasta ${n} ${seed}
    """
}

process runLCUBE {
    publishDir "$workflow.projectDir/../sampled/", mode: 'copy'
    // Note: long-term this needs to be updated to include a docker image for reproducibility purposes
    //        for now, set up your own R environment with the necessary packages.

    input:
    path seqs
    val seed
    val n
    val prefix
    path metadata

    output:
    path "${prefix}-${seed}-lcube.fasta", emit: seq
    path "${prefix}-${seed}-lcube.csv", emit: meta

    script:
    """
    Rscript $workflow.projectDir/scripts/lcube.r ${metadata} id Collection_Date ${seqs} ${prefix}-${seed}-lcube.csv ${prefix}-${seed}-lcube.fasta ${n} ${seed}
    """
}

process generateLphyScripts {
    publishDir "$workflow.projectDir/../lphy/", mode: 'copy'
    input:
    path seqs
    val seed
    val prefix

    output:
    path "${prefix}-${seed}.lphy", emit: lphy_script

    script:
    """
    sed "s|DATAGOESHERE|${subsample_fasta}|g" $workflow.projectDir/scripts/discrete_symmetric.lphy > ${prefix}-${seed}.lphy
    """
}

// process runBEAST {
//     // Note: long-term this needs to be updated for a BEAST 2 docker image with the requisite lphy package, but for now, set up your own BEAST 2 environment.

//     // This process uses LPHY to generate an xml file for BEAST 2, then runs BEAST 2 on the xml file.
//     input:
//     path seqs
//     val seed
//     val prefix

//     output:
//     path "${prefix}-${seed}-beast.xml", emit: beast_xml
//     path "${prefix}-${seed}-beast.log", emit: beast_log
//     path "${prefix}-${seed}-beast.trees", emit: beast_trees

//     script:
//     """
//     lphy -i ${seqs} -o ${prefix}-${seed}-beast.xml -t 1 -s ${seed}
//     beast -beagle -overwrite ${prefix}-${seed}-beast.xml > ${prefix}-${seed}-beast.log
//     """
// }

workflow{
    // Define the workflow for the simulation.

    // Define the parameters for the simulation.
    taxa = 1000
    seed = 12
    prefix = "sim_20240529"

    // Run the processes for the simulation.
    trueTree(taxa, seed, prefix)
    simulateSequences(trueTree.out.ttree, seed)
    
    // Simulate sequencing event, biased and standard.
    p=[0.8, 0.2, 0.8, 0.5, 0.2] // for now, using a hard-coded set of proportions. 
    s = runSimpleSampling(simulateSequences.out, seed, 500, prefix, trueTree.out.tmeta)
    simulateBiasedSampling(simulateSequences.out, trueTree.out.tmeta, seed, prefix, p)

    // Subsequence dataset using SRS, stratified, and LCUBE
    n = 100

    // Run subsampling on simple sequencing dataset
    simple_l = runLCUBE(s.seq, seed, n, prefix+"_simple", s.meta)
    simple_s = runSimpleSampling(s.seq, seed, n, prefix+"_simple", s.meta)
    simple_st = runStratifiedSampling(s.seq, seed, n, prefix+"_simple", s.meta)

    // Run subsampling on biased sequencing dataset
    biased_l = runLCUBE(simulateBiasedSampling.out.biased, seed, n, prefix+"_biased", simulateBiasedSampling.out.biased_meta)
    biased_s = runSimpleSampling(simulateBiasedSampling.out.biased, seed, n, prefix+"_biased", simulateBiasedSampling.out.biased_meta)
    biased_st = runStratifiedSampling(simulateBiasedSampling.out.biased, seed, n, prefix+"_biased", simulateBiasedSampling.out.biased_meta)

    // Generate LPHY scripts for BEAST 2
    generateLphyScripts(simple_l.out.seq, seed, prefix+"_simple_lcube")
    generateLphyScripts(simple_s.out.seq, seed, prefix+"_simple_simple")
    generateLphyScripts(simple_st.out.seq, seed, prefix+"_simple_strat")
    generateLphyScripts(biased_l.out.seq, seed, prefix+"_biased_lcube")
    generateLphyScripts(biased_s.out.seq, seed, prefix+"_biased_simple")
    generateLphyScripts(biased_st.out.seq, seed, prefix+"_biased_strat")

    // Run BEAST 2 on the subsampled datasets
}
