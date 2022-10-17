#!/bin/bash
#SBATCH --export=NONE
#SBATCH --mail-type=END
#SBATCH --mail-user=kdunton@uri.edu


# usage: andromeda
#        sbatch /home/kdunton/git/memonet/scripts/checksum.sh

module load scCATCH/3.1-foss-2021b-R-4.2.0 

Rscript /home/kdunton/git/memonet/scripts/checksum.r
