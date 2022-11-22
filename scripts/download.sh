#!/bin/bash

# download sequencing data from each mouse

wget --no-check-certificate -nH -np -r --cut-dirs=2 http://epigenomics.sdsc.edu/zhc268/lims_singlecell/slPsiwmg_JB_262_1_2_3/
wget --no-check-certificate -nH -np -r --cut-dirs=2 http://epigenomics.sdsc.edu/zhc268/lims_singlecell/QF9GKKgY_JB_263_1_2_3/
wget --no-check-certificate -nH -np -r --cut-dirs=2 http://epigenomics.sdsc.edu/jbuchanan/Komiyama/cellranger/JB_276_1_2_3/outs/
wget --no-check-certificate -nH -np -r --cut-dirs=2 http://epigenomics.sdsc.edu/jbuchanan/Komiyama/cellranger/JB_277_1_2_3/outs/
wget --no-check-certificate -nH -np -r --cut-dirs=2 http://epigenomics.sdsc.edu/jbuchanan/Komiyama/cellranger/JB_278_1_2_3/outs/
wget --no-check-certificate -nH -np -r --cut-dirs=2 http://epigenomics.sdsc.edu/zhc268/lims_singlecell/FHuOVucK_JB_279_1_2/
