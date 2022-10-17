#!/bin/bash
#SBATCH --export=NONE
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu


# usage: andromeda
#        sbatch /home/kdunton/git/memonet/scripts/download_fastqs.sh

wget --no-check-certificate -nH -np -r https://epigenomics.sdsc.edu/jbuchanan/Komiyama/fastqs/
  #--no-check-certificate: don't validate the site certificate
  #-np: don't ascend to parent directory
  #-nH don't make host directory (no idea what that means but Jing used it for the cellranger download)
