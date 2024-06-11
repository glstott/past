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

include {trueTree; simulateSequences; simulateBiasedSampling; runSimpleSampling; runSimpleSampling as simple2; runSimpleSampling as simple3; runLCUBE; runLCUBE as runLCUBE2; runStratifiedSampling; runStratifiedSampling as runStratifiedSampling2; generateLphyScripts as lphyLCUBE; generateLphyScripts as lphySimple; generateLphyScripts as lphyStratified;  generateLphyScripts as lphyLCUBE2; generateLphyScripts as lphySimple2; generateLphyScripts as lphyStratified2} from './simulation.nf'



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
//    simulateBiasedSampling(simulateSequences.out, seed, prefix, p, trueTree.out.tmeta)

    // Subsequence dataset using SRS, stratified, and LCUBE
    n = 100

    // Run subsampling on simple sequencing dataset
    simple_l = runLCUBE(s.seq, seed, n, prefix+"_simple", s.meta)
    simple_s = simple2(s.seq, seed, n, prefix+"_simple", s.meta)
    simple_st = runStratifiedSampling(s.seq, seed, n, prefix+"_simple", s.meta)

    // Run subsampling on biased sequencing dataset
 //   biased_l = runLCUBE2(simulateBiasedSampling.out.biased, seed, n, prefix+"_biased", simulateBiasedSampling.out.biased_meta)
 //   biased_s = simple3(simulateBiasedSampling.out.biased, seed, n, prefix+"_biased", simulateBiasedSampling.out.biased_meta)
 //   biased_st = runStratifiedSampling2(simulateBiasedSampling.out.biased, seed, n, prefix+"_biased", simulateBiasedSampling.out.biased_meta)

    // Generate LPHY scripts for BEAST 2
    lphyLCUBE(simple_l.seq, seed, prefix+"_simple_lcube")
    lphySimple(simple_s.seq, seed, prefix+"_simple_simple")
    lphyStratified(simple_st.seq, seed, prefix+"_simple_strat")
 //   lphyLCUBE2(biased_l.seq, seed, prefix+"_biased_lcube")
 //   lphySimple2(biased_s.seq, seed, prefix+"_biased_simple")
 //   lphyStratified2(biased_st.seq, seed, prefix+"_biased_strat")

    // Run BEAST 2 on the subsampled datasets
}
