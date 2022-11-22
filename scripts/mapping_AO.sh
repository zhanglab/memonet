#!/bin/bash
#SBATCH --mem=250G 
##SBATCH -c 4
#SBATCH --export=NONE
##SBATCH --exclusive



# usage: 
#        sbatch <PATH>/mapping_AO.sh
# takes ~7min


# load Seurat and other dependencies indicated in R script



Rscript <PATH>/mapping_AO.r
