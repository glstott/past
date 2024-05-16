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

process simulateSequences {
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
    //simulate the acquisition of a convenience sample from the simulated sequences
    input:
    path seqs
    val seed
    val prefix
    val p

    output:
    path "${prefix}-${seed}-convenience.fa", emit: simple

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

process runLCUBE {
    // Note: long-term this needs to be updated to include a docker image for reproducibility purposes
    //        for now, set up your own R environment with the necessary packages.

    input:
    path seqs
    val seed
    val n
    val prefix
    path metadata

    output:
    path "${prefix}-${seed}-lcube.fasta", emit: lcube
    path "${prefix}-${seed}-lcube.csv", emit: lcube_meta

    script:
    """
    Rscript $workflow.projectDir/scripts/lcube.r ${metadata} id Collection_Date ${seqs} ${prefix}-${seed}-lcube.csv ${prefix}-${seed}-lcube.fasta ${n} ${seed}
    """
}

process generateLphyScripts {
    // Note: long-term this needs to be updated to include a docker image for reproducibility purposes
    //        for now, set up your own R environment with the necessary packages.

    input:
    path seqs
    val seed
    val prefix

    output:
    path "${prefix}-${seed}.lphy", emit: lphy_script

    shell:
    """
    cat > !{prefix}-!{seed}.lphy << EOL
    data {
    // Specify options for reading the date information from file
    options = {ageDirection="forward", ageRegex=".*_.*_(\d*\.\d+|\d+\.\d*)$"};
    D = readNexus(file="!{seqs}", options=options);

    // Retrieve number of sites in alignment and taxa names
    L = D.nchar();
    taxa = D.taxa();

    // Extract Trait information and count # of traits
    D_trait = extractTrait(taxa=taxa, sep="_", i=0);
    K = D_trait.canonicalStateCount();

    // Specify the data type for traits and number of pairwise combos
    dim = K*(K-1)/2;
    dataType = D_trait.dataType();
    }
    model {
    // Model evolutionary history
    Π ~ Dirichlet(conc=[2.0, 2.0, 2.0, 2.0]);      // Nucleotide Frequency prior
    κ ~ LogNormal(meanlog=1.0, sdlog=1.25);        // Transition/transversion ratio prior
    γ ~ LogNormal(meanlog=0.0, sdlog=2.0);         // Shape parameter for discretized gamma
    r ~ DiscretizeGamma(shape=γ, ncat=4, 
            replicates=L);                           // Site rates
    Θ ~ LogNormal(meanlog=0.0, sdlog=1.0);         // Effective population size prior
    ψ ~ Coalescent(taxa=taxa, theta=Θ);            // Coalescent time scaled phylogenetic tree prior
    D ~ PhyloCTMC(Q=hky(kappa=κ, freq=Π), mu=0.004, 
            siteRates=r, tree=ψ);                    // MCMC for evolution of alignment

    // Model the geographic history
    // modeled via P(t) = e^(Dt)
    π_trait ~ Dirichlet(conc=rep(element=3.0, times=K));   // Base frequencies, assuming symmetry
    I ~ Bernoulli(minSuccesses=dim-2, p=0.5, 
            replicates=dim);                                 // Determines which rates are 0. 
                                                            // This, plus select, implements BSSVS.
    R_trait ~ Dirichlet(conc=rep(element=1.0, times=dim)); // Off-diagonal entries of the rate matrix.
    Q_trait = generalTimeReversible(rates=select(x=R_trait, 
                indicator=I), freq=π_trait);               // intantaneous rate matrix
    μ_trait ~ LogNormal(meanlog=0, sdlog=1.25);            // migration events per unit time
    D_trait ~ PhyloCTMC(L=1, Q=Q_trait, dataType=dataType, 
                mu=μ_trait, tree=ψ);                       // MCMC for discrete locations
    }
    EOL
    """
}

process runBEAST {
    // Note: long-term this needs to be updated for a BEAST 2 docker image with the requisite lphy package, but for now, set up your own BEAST 2 environment.

    // This process uses LPHY to generate an xml file for BEAST 2, then runs BEAST 2 on the xml file.
    input:
    path seqs
    val seed
    val prefix

    output:
    path "${prefix}-${seed}-beast.xml", emit: beast_xml
    path "${prefix}-${seed}-beast.log", emit: beast_log
    path "${prefix}-${seed}-beast.trees", emit: beast_trees

    script:
    """
    lphy -i ${seqs} -o ${prefix}-${seed}-beast.xml -t 1 -s ${seed}
    beast -beagle -overwrite ${prefix}-${seed}-beast.xml > ${prefix}-${seed}-beast.log
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

    // Simulate sequencing event, biased and standard.
    p=[0.8, 0.2, 0.8, 0.5, 0.1] // for now, using a hard-coded set of proportions. 
    simulateSimpleSampling(simulateSequences.out, seed, prefix, 0.1)
    simulateBiasedSampling(simulateSequences.out, seed, prefix, p)

    // Subsequence dataset using SRS, stratified, and LCUBE
}
