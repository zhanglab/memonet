This file details all steps of the analysis for Chapter 1, including scripts and location of input/output data files.

# Set-up
## Set up working directory structure
```{r} 
cd ~/Downloads
mkdir RNAseq
  # this will be the overarching directory
cd RNAseq
mkdir data
cd data
mkdir our_data AIBS_data
```

## Download Packages 
Cell Ranger
- Follow steps here (make note of download location): https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

Seurat


DESC


## Download datasets
Our dataset:
- Location of data: https://docs.google.com/spreadsheets/d/1mU7l8Oj-Fr4FYE6IlmTcX_s2YdCtZXx0GgJI3jBE3J4/edit#gid=0 (* will need to change this to the ncbi repository number once we upload data there)
- Script for downloading: /data/zhanglab/jingwang/brain/RNAseq/Takaki/Deep/download.sh
- Download to: RNAseq/data/our_data/  

AIBS dataset:
- Location of data:
- Download to: RNAseq/data/AIBS_data

# Pre-processing
## Aggregate the data of all mice
Use the aggr command to combine the data of all mice
- Wd: RNAseq/data/our_data
- Generate a libraries.csv file with the locations of each mouse's data, in the following format: 
![](embedded_images/libraries.csv.png)
  - the order of rows determines the suffix attached to the barcodes of each mouse (in order to differentiate which mouse is which in the combined file), ie mouse 262 will have -1, 263 -2, etc
  - trained mice: 263,277,279 ie -2,-4,-6
  - control mice: 262,276,278 ie -1,-3,-5
- Run aggr command:
  - First run this line on the terminal so that Cell Ranger package can be found, replacing the path with the location of your Cell Ranger download location: 
  ```{r} 
  export PATH=<path to cellranger download>/cellranger-4.0.0:$PATH
  ```
  - Then run aggr. --none turns off depth normalization, due to the requirement of the DESC clustering package needing unnormalized counts as input. --id is the name of the output file
  ```{r}
  nohup srun -o jobs%j.out -c 20 \cellranger aggr --csv=/data/our_data/libraries.csv --none --id=combined_cellranger_no-normalization &
  ```
   - The output directory will be found here: RNAseq/data/our_data/combined_cellranger_no-normalization/outs/filtered_feature_bc_matrix
 
## Determine QC thresholds
```{r}
cd ~/Downloads/RNAseq
mkdir QC
```

Script: /data/zhanglab/kdunton/neuron_model/katie-scripts/RNAseq/scripts_unorganized/QCthresholds.r 
- Choose cutoffs based on the plots with no cutoff, then update code to use those cutoffs and visualize the remaining dataset

Usage: sbatch /data/zhanglab/kdunton/neuron_model/katie-scripts/RNAseq/scripts_unorganized/QCthresholds.sh

Wd: RNAseq/QC

Input: 
- Data directory for each mouse, ie /data/zhanglab/jingwang/brain/RNAseq/Takaki/Deep/deepseq_2/slPsiwmg_JB_262_1_2_3/filtered_feature_bc_matrix/

Output:
- ctrl_without_cutoff.png: violin plots of the 3 control mice before QC
- QC_ctrl.png: violin plots of the 3 control mice after QC
- train_without_cutoff.png
- QC_train.png


# Map our data onto AIBS dataset ("OA mapping")
```{r}
cd ~/Downloads/RNAseq
mkdir AIBSmapping
cd AIBSmapping
mkdir OA test
  # 'OA' refers to Our data mapped to AIBS data
```

Script: pairwiseOA_clean.r

Wd: RNAseq/AIBSmapping/OA

Query: our data
Reference: AIBS data

Input: 
- 10X directory of all mice combined: RNAseq/data/our_data/combined_cellranger_no-normalization/outs/filtered_feature_bc_matrix/
- AIBS dataset. The three files here: RNAseq/data/AIBS_data

Output: 
- prediction_scores.csv: lists each cell, the predicted cell type label, and prediction scores for each cell type
- umap_referenceOA-celltypes.png: umap of the reference (AIBS) with cell type labels
- umap_OA_query-predictedLabels.png: umap of the query (our data) projected onto AIBS space, labeled with the predicted labels

## Test the accuracy of mapping on AIBS data
```{r}
cd ~/Downloads/RNAseq/AIBSmapping/test
```

### 1. Generate a test dataset (downsample to 25% of each cell type; remove sample cells from rest of reference) and perform label transfer from the remaining 75% of data. Do this 100 times.

Query: 25% of AIBS data

Reference: remaining 75% of AIBS data

Script: testA_prediction_cutoff.r

Wd: RNAseq/AIBSmapping/test

Input: 
- AIBS dataset. The three files here: RNAseq/data/AIBS_data
 
Output:
- 100 prediction_scores_*.csv files

### 2. Combine the 100 prediction score files into one file
Script: testA_prediction_cutoff_table.r

Wd: RNAseq/AIBSmapping/test

Input: 
- 100 prediction_scores_*.csv files
 
Output:
- prediction_cutoff.csv

### 3. Summarize results in RStudio: calculate false classification percentage, generate confusion matrix, compare mean and median scores of the test with OA mapping
Script: testA_prediction_cutoff2.r

Wd: RNAseq/AIBSmapping/test

Input: 
- prediction_cutoff.csv
- ~/Downloads/RNAseq/AIBSmapping/OA/prediction_scores.csv

Output:
- AIBStest_confusionMtx.png: heatmap of original cell type labels vs predicted labels for AIBS testing
- maxPredictionScores-AIBStest_and_OA.csv: table of mean and median prediction.score.max for AIBS testing and OA mapping; shows calculation for all celltypes and subset for L2/3 cells

# Investigate cell type proportions

script: OverallDatasetDescription2.r

functions:
- calculate cell subclass proportion
- calculate inhibitory neuron percentage
- calculate neuron subclass proportion (neurons only)
- compare neuron subclass proportion in our dataset and AIBS
- 
input:
- AIBS metadata: ~/Downloads/RNAseq/data/AIBS_data/aibs_barcodes.tsv
- cell type predictions for our data: ~/Downloads/RNAseq/AIBSmapping/OA/prediction_scores.csv



