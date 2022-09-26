#!/bin/bash
#SBATCH --export=NONE
#SBATCH -p uri-cpu,cpu,cpu-preempt


#usage: sbatch --mem=120G /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESeq2-table.sh


module load miniconda
conda activate DESeq2


Rscript /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESeq2-table.r
