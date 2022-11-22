#!/bin/bash
#SBATCH --export=NONE

#it won't run correctly with multiple nodes if you're using BiocParallel


#usage: 
#       sbatch --exclusive --mem=250G <PATH>/DESeq2_tr_vs_ctrl.sh
#       takes 16 min

# load miniconda
# load DESeq2 and other dependencies listed in R script


export OMP_NUM_THREADS=1

Rscript <PATH>/DESeq2_tr_vs_ctrl.r
