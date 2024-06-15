#!/bin/bash
#SBATCH --job-name=beast-batches
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=%j.out

ml Beast/2.7.5-foss-2022a-CUDA-11.7.0

beast -beagle_info
beast -beagle -seed 12 WHICH
