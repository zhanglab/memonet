#!/bin/bash
#SBATCH -p uri-cpu,cpu,cpu-preempt
#SBATCH -c 2
#SBATCH --export=NONE
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu


#it won't run correctly with multiple nodes if you're using BiocParallel

#usage: ANDROMEDA or bluewaves
#       sbatch --exclusive --mem=120G /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESeq2_clusterX_vs_clusterY_whole.sh
#         takes ~20min  on andromeda and unity

#module load R-bundle-Bioconductor/3.10-foss-2019b
#module load scCATCH/2.1-foss-2019b-R-3.6.2
module load miniconda
conda activate DESeq2

export OMP_NUM_THREADS=1

Rscript /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESeq2_clusterX_vs_clusterY_whole.r
