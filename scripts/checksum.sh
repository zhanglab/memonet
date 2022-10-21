#!/bin/bash
#SBATCH --export=NONE
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu


# usage: andromeda
#        sbatch /home/kdunton/git/memonet/scripts/checksum.sh

module load scCATCH/2.1-foss-2020b-R-4.0.3 

Rscript /home/kdunton/git/memonet/scripts/checksum.r
