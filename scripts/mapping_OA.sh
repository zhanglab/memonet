#!/bin/sh
#SBATCH --mem=250G 
##SBATCH -c 4
#SBATCH --export=NONE
#SBATCH --exclusive


# usage: 
#        sbatch <PATH>/mapping_OA.sh
# takes ~ 1hr at 120G, 30min at 250G


# load Seurat and other dependencies listed in R script

Rscript <PATH>/mapping_OA.r
