#!/bin/bash
#SBATCH --mem=250G
#SBATCH --export=NONE

# usage: sbatch ~/Downloads/RNAseq/memonet_github_repo/memonet/scripts/QCthresholds.sh

module load scCATCH/2.1-foss-2019b-R-3.6.2
  #andromeda


Rscript ~/Downloads/RNAseq/memonet_github_repo/memonet/scripts/QCthresholds.r
