#!/bin/bash
#SBATCH --mem=120G
#SBATCH -p uri-cpu,cpu,cpu-preempt
#SBATCH --export=NONE
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu

# using -q zhanglab means it can't use n107

# usage: andromeda
#        sbatch /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/labelTransfer_AO_umapDESC.sh
# takes ~7min

#module load scCATCH/2.1-foss-2019b-R-3.6.2
#module load scCATCH/2.1-foss-2020b-R-4.0.3
  #this is Seurat v4
module load uri
module load scCATCH/3.1-foss-2021b-R-4.2.0
  #this will automatically load R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0

Rscript /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/labelTransfer_AO_umapDESC.r
