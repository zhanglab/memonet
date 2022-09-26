#!/bin/bash
#SBATCH --export=NONE
#SBATCH -c 2
#SBATCH -p uri-cpu,cpu
#SBATCH --mem=250G
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu

#usage: sbatch /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/AvgGeneExp.sh
# takes 16hrs for all genes, was only using 1 cpu I think

module load miniconda
conda activate DESeq2

Rscript /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/AvgGeneExp.r
