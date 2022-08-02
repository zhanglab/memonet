This file details all steps of the analysis, including scripts and location of input/output data files.

# Raw snRNA-seq data processing
Raw reads were processed using the Cell Ranger v4.0.0 pipeline.

The following steps were performed by scientists at 10X: 
- the mkfastq command was used to demultiplex the BCL files and create fastq files
- the count command was used, with parameter --expect-cells set to 7000, to perform alignment, barcode and UMI counting to generate feature-barcode matrices

Feature barcode matrices were downloaded:
- Location of feature barcode matrices: https://docs.google.com/spreadsheets/d/1mU7l8Oj-Fr4FYE6IlmTcX_s2YdCtZXx0GgJI3jBE3J4/edit#gid=0 
- Script for downloading: /data/zhanglab/jingwang/brain/RNAseq/Takaki/Deep/download.sh
- Download locations:  
  - 262,263,279 had enough data at the second round and can be found here: /data/zhanglab/jingwang/brain/RNAseq/Takaki/Deep/deepseq_2
  - 276,277,278 were sequenced a third round and can be found here: /data/zhanglab/jingwang/brain/RNAseq/Takaki/Deep/deepseq_3

The aggr command was used to aggregate the data for all mice. The parameter --none was used to turn off depth normalization, due to the requirement of the DESC clustering package needing unnormalized counts as input.
