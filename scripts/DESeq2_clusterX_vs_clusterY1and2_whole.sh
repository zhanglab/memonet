#!/bin/bash
#SBATCH -p uri-cpu,cpu,cpu-preempt
##SBATCH -t 100:00:00
##SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
##SBATCH -q zhanglab
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu


#zhanglab partition moves you to high mem node even without setting mem=500G
#-t says to give maximum 100 hours
#This specifies to allocate 1 node and 48 processors (cores) per node
  #set to all 48 cores so you take up the whole node (same as using --exclusive)
#it won't run correctly with multiple nodes if you're using BiocParallel

#usage: ANDROMEDA or bluewaves
#       sbatch --exclusive --mem=120G /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESeq2_clusterX_vs_clusterY1and2_whole.sh
#         takes ~20min  on andromeda and unity

#module load R-bundle-Bioconductor/3.10-foss-2019b
#module load scCATCH/2.1-foss-2019b-R-3.6.2
module load miniconda
conda activate DESeq2

export OMP_NUM_THREADS=1

Rscript /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESeq2_clusterX_vs_clusterY1and2_whole.r
