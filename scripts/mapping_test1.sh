#!/bin/bash
#SBATCH --mem=30G
  #each job takes this much memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9
  #this means 9 threads
#SBATCH --array=1-100%20
#SBATCH --export=NONE



# usage: 
#        sbatch  <PATH>/mapping_test1.sh


# load Seurat and other dependencies listed in R script


Rscript  <PATH>/mapping_test1.r $SLURM_ARRAY_TASK_ID
