#!/bin/bash
#SBATCH --mem=250G
#SBATCH -p uri-cpu,cpu,cpu-preempt 
#SBATCH --export=NONE
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu


# usage: andromeda
#        sbatch /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/pairwiseOA_clean.sh
# takes ~ 1hr at 120G, 30min at 250G


source /modules/lmod/lmod/lmod/init/sh
module load uri 
module load scCATCH/3.1-foss-2021b-R-4.2.0
  #this will automatically load R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0

Rscript /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/pairwiseOA_clean.r
