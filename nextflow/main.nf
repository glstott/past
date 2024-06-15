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

params.n = 500
params.taxa = 5000
params.sequenced = 2500
params.prefix = "sim_0613"
params.p = "c(0.1,0.8,0.8,0.8,0.8)"
params.loc="4"


process trueTree {
    //conda "-c bioconda -c conda-forge r-base r-ape=5.7 r-diversitree r-phangorn=2.11 bioconductor-treeio=1.26.0"
    //conda "$workflow.projectDir/envs/simulation.yml"
    conda "-c bioconda iqtree"
    publishDir "$workflow.projectDir/../raw/", mode: 'copy'
    // Process for generating true trees given a number of taxa.

    input:
    val seed
    val taxa
    val prefix

    output:
    tuple path("${prefix}-${seed}.fa"), path("${prefix}-${seed}.csv"), val(seed), emit: data
    path "${prefix}-${seed}.nwk", emit: ttree
    path "${prefix}-${seed}.nex", emit: ttree_nex


    script:
    """
    Rscript $workflow.projectDir/scripts/simulate.r --vanilla ${prefix}-${seed} ${prefix}-${seed}.csv ${taxa} ${seed}
    iqtree2 --alisim ${prefix}-${seed} -m HKY --branch-scale 0.004 --length 1500 -t ${prefix}-${seed}.nwk --out-format fasta  ${seed}
    """
}

process convenienceSampling {
    publishDir "$workflow.projectDir/../raw/", mode: 'copy'
    input:
    tuple path(seqs), path(metadata), val(seed)
    val n

    output:
    tuple path("${seqs.simpleName}_s${n}.fasta"), path("${seqs.simpleName}_s${n}.csv"), val(seed), emit: data

    script:
    """
    Rscript $workflow.projectDir/scripts/convenienceSample.R --vanilla ${metadata} id Collection_Date ${seqs} ${seqs.simpleName}_s${n}.csv ${seqs.simpleName}_s${n}.fasta ${n} ${seed}
    """
}

process geoBiasedSampling {
    publishDir "$workflow.projectDir/../raw/", mode: 'copy'
    input:
    tuple path(seqs), path(metadata), val(seed)
    val n
    val p

    output:
    tuple path("${seqs.simpleName}_b${n}.fasta"), path("${seqs.simpleName}_b${n}.csv"), val(seed), emit: data

    script:
    """
    Rscript $workflow.projectDir/scripts/biasedSample.R --vanilla ${metadata} id Collection_Date ${seqs} ${seqs.simpleName}_b${n}.csv ${seqs.simpleName}_b${n}.fasta ${n} ${seed} ${p}
    """
}

process temporalBiasedSampling {
    publishDir "$workflow.projectDir/../raw/", mode: 'copy'
    input:
    tuple path(seqs), path(metadata), val(seed)
    val n
    val loc

    output:
    tuple path("${seqs.simpleName}_${loc}t${n}.fasta"), path("${seqs.simpleName}_${loc}t${n}.csv"), val(seed), emit: data

    script:
    """
    Rscript $workflow.projectDir/scripts/temporalBiasedSample.R --vanilla ${metadata} id Collection_Date ${seqs} ${seqs.simpleName}_${loc}t${n}.csv ${seqs.simpleName}_${loc}t${n}.fasta ${n} ${seed} ${loc}
    """
}

process runSimpleSampling {
    publishDir "$workflow.projectDir/../lphy/sampled/", mode: 'copy'
    input:
    tuple path(seqs), path(metadata), val(seed)
    val n

    output:
    tuple path("${seqs.simpleName}_s${n}.fasta"), path("${seqs.simpleName}_s${n}.csv"), val(seed), emit: data

    script:
    """
    Rscript $workflow.projectDir/scripts/simpleSample.R --vanilla ${metadata} id Collection_Date ${seqs} ${seqs.simpleName}_s${n}.csv ${seqs.simpleName}_s${n}.fasta ${n} ${seed}
    """
}

process runStratifiedSampling {
    publishDir "$workflow.projectDir/../lphy/sampled/", mode: 'copy'
    input:
    tuple path(seqs), path(metadata), val(seed)
    val n

    output:
    tuple path("${seqs.simpleName}_t${n}.fasta"), path("${seqs.simpleName}_t${n}.csv"), val(seed), emit: data

    script:
    """
    Rscript $workflow.projectDir/scripts/stratifiedSample.R --vanilla ${metadata} id Collection_Date ${seqs} ${seqs.simpleName}_t${n}.csv ${seqs.simpleName}_t${n}.fasta ${n} ${seed}
    """
}

process runLCUBE {
    publishDir "$workflow.projectDir/../lphy/sampled/", mode: 'copy'
    // Note: long-term this needs to be updated to include a docker image for reproducibility purposes
    //        for now, set up your own R environment with the necessary packages.
    input:
    tuple path(seqs), path(metadata), val(seed)
    val n

    output:
    tuple path("${seqs.simpleName}_l${n}.fasta"), path("${seqs.simpleName}_l${n}.csv"), val(seed), emit: data

    script:
    """
    Rscript $workflow.projectDir/scripts/lcube.r --vanilla ${metadata} id Collection_Date ${seqs} ${seqs.simpleName}_l${n}.csv ${seqs.simpleName}_l${n}.fasta ${n} ${seed}
    """
}

process generateLphyScripts {
    publishDir "$workflow.projectDir/../lphy/", mode: 'copy'
    input:
    tuple path(seqs), path(metadata), val(seed)

    output:
    path "${seqs.simpleName}.lphy", emit: lphy_script
    path "${seqs.simpleName}-100M.xml", emit: beast_xml

    script:
    """
    sed "s|DATAGOESHERE|${seqs}|g" $workflow.projectDir/scripts/discrete_symmetric.lphy > ${seqs.simpleName}.lphy
    lphybeast ${seqs.simpleName}.lphy -l 100000000 -o ${seqs.simpleName}-100M.xml
    """
}


workflow {
    c=channel.of(1..5) 
    t=trueTree(seed=c, taxa=params.taxa, prefix=params.prefix)
    simpleConvenience=convenienceSampling(t.data, params.sequenced)
    biasedConvenience=geoBiasedSampling(t.data, params.sequenced, params.p)
    temporalBiasedConvenience=temporalBiasedSampling(t.data, params.sequenced, params.loc)
    c=simpleConvenience.mix(biasedConvenience,temporalBiasedConvenience)
    simple=runSimpleSampling(c.data, params.n)
    stratified=runStratifiedSampling(c.data, params.n)
    lcube=runLCUBE(c.data, params.n)
    simple.mix(stratified,lcube) | generateLphyScripts
}
