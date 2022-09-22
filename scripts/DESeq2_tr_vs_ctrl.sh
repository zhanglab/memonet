#!/bin/bash
#SBATCH -p uri-cpu,cpu,cpu-preempt
#SBATCH --export=NONE
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu

#it won't run correctly with multiple nodes if you're using BiocParallel


#usage: ANDROMEDA or bluewaves (might work at 120G)
#       sbatch --exclusive --mem=250G /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESeq2_tr_vs_ctrl.sh
#       takes 16 min on unity with 250G

module load miniconda
conda activate DESeq2


export OMP_NUM_THREADS=1

Rscript /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESeq2_tr_vs_ctrl.r
