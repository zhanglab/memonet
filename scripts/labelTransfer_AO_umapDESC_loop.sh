#!/bin/bash
#SBATCH --mem=250G 
##SBATCH -c 4
#SBATCH --export=NONE
#SBATCH -q zhanglab
##SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu

# using -q zhanglab means it can't use n107

# usage: andromeda
#        sbatch /home/kdunton/git/memonet/scripts/labelTransfer_AO_umapDESC_loop.sh
# takes ~7min

#module load scCATCH/2.1-foss-2019b-R-3.6.2
module load scCATCH/2.1-foss-2020b-R-4.0.3
  #this is Seurat v4

#Rscript /home/kdunton/git/memonet/scripts/labelTransfer_AO_umapDESC_loop.r
Rscript /data/zhanglab/kdunton/neuron_model/katie-scripts/RNAseq/BICCNscripts/labelTransfer_AO_umapDESC_loop.r
