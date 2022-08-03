This file details all steps of the analysis for Chapter 1, including scripts and location of input/output data files.

# Set-up
## Set up working directory structure
```{r} 
mkdir RNAseq
cd RNAseq
mkdir data QC
cd data
mkdir our_data AIBSmapping
```

## Download Packages 
Cell Ranger
- Follow steps here (make note of download location): https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

Seurat


DESC


## Download datasets
Our dataset:
- Location of data: https://docs.google.com/spreadsheets/d/1mU7l8Oj-Fr4FYE6IlmTcX_s2YdCtZXx0GgJI3jBE3J4/edit#gid=0 
- Script for downloading: /data/zhanglab/jingwang/brain/RNAseq/Takaki/Deep/download.sh
- Download to: data/our_data/  

AIBS dataset:
- Location of data:
- Download to: data/AIBS_data

# Pre-processing
## Aggregate the data of all mice
Use the aggr command to combine the data of all mice
- Wd: /data/our_data
- Generate a libraries.csv file with the locations of each mouse's data, in the following format: 
![](embedded_images/libraries.csv.png)
  - the order of rows determines the suffix attached to the barcodes of each mouse (in order to differentiate which mouse is which in the combined file), ie mouse 262 will have -1, 263 -2, etc
- Run aggr command:
  - First run this line on the terminal so that Cell Ranger package can be found, replacing the path with the location of your Cell Ranger download location: 
  ```{r} 
  export PATH=<path to cellranger download>/cellranger-4.0.0:$PATH
  ```
  - Then run aggr. --none turns off depth normalization, due to the requirement of the DESC clustering package needing unnormalized counts as input. --id is the name of the output file
  ```{r}
  nohup srun -o jobs%j.out -c 20 \cellranger aggr --csv=/data/our_data/libraries.csv --none --id=combined_cellranger_no-normalization &
  ```
   - The output directory will be found here: /data/our_data/combined_cellranger_no-normalization/outs/filtered_feature_bc_matrix
 




