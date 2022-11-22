#!/bin/bash
#SBATCH --gpus=2
#SBATCH -c 2 
#SBATCH --mem=250g


#--gpus parameter requests access to a GPU node with 2 GPUs
# -c parameter requests 2 CPU cores


#usage: 
#       sbatch <PATH>/DESC_EDGclustering.sh
# takes ~ 5min


# load libcudart.so.11.0
# load DESC and other dependencies listed in the Python script


python <PATH>/DESC_EDGclustering.py
