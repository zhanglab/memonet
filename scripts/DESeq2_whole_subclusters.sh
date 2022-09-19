#!/bin/bash
##SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu



#This specifies to allocate 1 node and 48 processors (cores) per node
  #set to all 48 cores so you take up the whole node (same as using --exclusive)
#it won't run correctly with multiple nodes if you're using BiocParallel

#usage: ANDROMEDA or bluewaves (might work at 120G)
#       sbatch --exclusive --mem=250G /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESeq2_whole_subclusters.sh
#         takes ~1hr 40min on andromeda

module load miniconda
conda activate DESeq2


export OMP_NUM_THREADS=1

Rscript /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESeq2_whole_subclusters.r
