#!/bin/bash
#SBATCH --gpus=2 #set to 1 when Kevin tested
#SBATCH -c 2 #set to 1 when Kevin tested
#SBATCH --mem=250g
#SBATCH --export=NONE
#SBATCH -p gpu,gpu-preempt
  # tells which type of nodes you want to use 
#SBATCH --time=4:00:00
##SBATCH --nodelist=
##SBATCH -q zhanglab
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu

#--gpus-per-node parameter requests access to a GPU node with 2 GPUs
# -c parameter requests 2 CPU cores


#usage: RUN ON ANDROMEDA, NOT BLUEWAVES
#       sbatch /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESC_DGclustering.sh
# takes ~ 5min

module load cuda
 # for libcudart.so.11.0
module load cudnn/cuda11-8.4.1.50
module load miniconda
  # to activate env
#conda activate desc_clone
conda activate /work/pi_yingzhang_uri_edu/.conda/envs/desc

#module load desc/2.1.1-fosscuda-2020b
  #there was an error ' 'str object has no attribute decode' someitmes when running script, and other times it would run fine without changing anything; appears to be an issue with tensorflow (issue should have been with 2.1 but this module uses 2.4); Kevin installed 2.5 in new module
#module load desc/2.1.1-foss-2021a
  #started using 1/21/22

python /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESC_DGclustering.py
