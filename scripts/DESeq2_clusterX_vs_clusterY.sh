#!/bin/bash
#SBATCH -c 2
#SBATCH --export=NONE


#it won't run correctly with multiple nodes if you're using BiocParallel

#usage:
#       sbatch --exclusive --mem=120G <PATH>/DESeq2_clusterX_vs_clusterY.sh
#         takes ~20min


# load miniconda
# load DESeq2 and other dependencies listed in the R script

export OMP_NUM_THREADS=1

Rscript <PATH>/DESeq2_clusterX_vs_clusterY.r
