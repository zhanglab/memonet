#!/bin/bash
#SBATCH --mem=250G
#SBATCH --export=NONE

# usage: sbatch <PATH>/QCthresholds.sh


# load Seurat


Rscript <PATH>/QCthresholds.r
