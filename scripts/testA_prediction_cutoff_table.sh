#!/bin/bash
#SBATCH --mem=250G
#SBATCH --export=NONE

# usage: andromeda or bluewaves
#        sbatch ~/Downloads/RNAseq/memonet_github_repo/memonet/scripts/testA_prediction_cutoff_table.sh
# takes < 5min


module load R/3.6.2-foss-2019b

Rscript ~/Downloads/RNAseq/memonet_github_repo/memonet/scripts/testA_prediction_cutoff_table.r
