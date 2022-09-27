#!/bin/bash
#SBATCH --gpus=2 
#SBATCH -c 2 
#SBATCH --mem=250g
#SBATCH --export=NONE
#SBATCH -p gpu,gpu-preempt
  # tells which type of nodes you want to use 
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu

#--gpus-per-node parameter requests access to a GPU node with 2 GPUs
# -c parameter requests 2 CPU cores


#usage: RUN ON ANDROMEDA, NOT BLUEWAVES
#       sbatch /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESC_highRankGenes_step1_loop15.sh
# takes ~ 5min

module load cuda
 # for libcudart.so.11.0
module load cudnn/cuda11-8.4.1.50
module load miniconda
  # to activate env
conda activate /work/pi_yingzhang_uri_edu/.conda/envs/desc


python /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESC_highRankGenes_step1_loop15.py
