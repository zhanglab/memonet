#!/bin/bash
#SBATCH --mem=30G
  #each job takes this much mem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9
  #this means 9 threads
#SBATCH --array=1-100%20
##SBATCH --array=1-20%20
  #this will put you on 5 nodes bc 20 jobs with 4 per node (30G each job, 128G each node)
#SBATCH --export=NONE
#SBATCH -q zhanglab
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu

# using -q zhanglab means it can't use n107

# usage: andromeda
#        sbatch  ~/Downloads/RNAseq/memonet_github_repo/memonet/scripts/testA-prediction_cutoff.sh

module load scCATCH/2.1-foss-2020b-R-4.0.3
  #this is Seurat v4

Rscript  ~/Downloads/RNAseq/memonet_github_repo/memonet/scripts/testA-prediction_cutoff.r $SLURM_ARRAY_TASK_ID
