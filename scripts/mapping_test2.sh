#!/bin/bash
#SBATCH --mem=250G
#SBATCH --export=NONE

# usage: 
#        sbatch <PATH>/mapping_test2.sh
# takes < 5min


# load R and dependencies listed in R script

Rscript <PATH>/mapping_test2.r
