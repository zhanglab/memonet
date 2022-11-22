#!/bin/bash
#SBATCH --export=NONE


#usage: sbatch --mem=120G <PATH>/DESeq2-table.sh


#load miniconda
#load DESeq2 and other dependences listed in the R script


Rscript <PATH>/DESeq2-table.r
