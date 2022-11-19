This file details all steps of the analysis for Chapter 1, including scripts and location of input/output data files.

# Set-up
## Set up working directory structure
```{r} 
cd ~/Downloads
mkdir RNAseq
  # this will be the project directory
cd RNAseq
mkdir data
cd data
mkdir memonet_data AIBS_data
```

## Clone repo
```{r} 
cd ~/Downloads/RNAseq
mkdir memonet_github_repo
cd memonet_github_repo
git clone git@github.com:zhanglab/memonet.git
```

You can find scripts here: ~/Downloads/RNAseq/memonet_github_repo/memonet/scripts


## Download Packages 
Cell Ranger
- Follow steps here (make note of download location): https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

Seurat *

DESC *

DESeq2 *

EnrichR *


## Download datasets
```{r} 
cd ~/Downloads/RNAseq/memonet_data
mkdir data_download
```

MEMONET dataset:
- Location of data: https://docs.google.com/spreadsheets/d/1mU7l8Oj-Fr4FYE6IlmTcX_s2YdCtZXx0GgJI3jBE3J4/edit#gid=0 (* will need to change this to the ncbi repository number once we upload data there)
- Script for downloading: /data/zhanglab/jingwang/brain/RNAseq/Takaki/Deep/download.sh
- Download to: ~/Downloads/RNAseq/data/memonet_data/data_download

AIBS dataset:
- Location of data:
- Download to: ~/Downloads/RNAseq/data/AIBS_data
- Rename the files to: aibs_barcodes.tsv, aibs_genes.tsv, aibs_matrix.mtx

# Pre-processing of memonet data
## Aggregate the data of all mice using Cell Ranger commands
Use Cell Ranger's aggr command to combine the data of all mice:
- Wd: ~/Downloads/RNAseq/data/memonet_data
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
  nohup srun -o jobs%j.out -c 20 \cellranger aggr --csv=/data/memonet_data/libraries.csv --none --id=combined_cellranger_no-normalization &
  ```
   - The output directory will be found here: ~/Downloads/RNAseq/data/memonet_data/combined_cellranger_no-normalization/outs/filtered_feature_bc_matrix
 
## Determine QC thresholds
```{r}
cd ~/Downloads/RNAseq
mkdir QC
```

Script: ~/Downloads/RNAseq/memonet_github_repo/memonet/scripts/QCthresholds.r 
- Thresholds were chosen based on violin plots of the data before QC

Wd: ~/Downloads/RNAseq/QC

Input: 
- Data directory for each mouse, ie ~/Downloads/RNAseq/data/memonet_data/data_download/slPsiwmg_JB_262_1_2_3/filtered_feature_bc_matrix/

Output:
- ctrl_without_cutoff.png: violin plots of the 3 control mice before QC
- QC_ctrl.png: violin plots of the 3 control mice after QC
- train_without_cutoff.png
- QC_train.png


# Map MEMONET data onto AIBS dataset ("OA mapping")
```{r}
cd ~/Downloads/RNAseq
mkdir AIBSmapping
cd AIBSmapping
mkdir OA test
  # 'OA' refers to Our (MEMONET) data mapped to AIBS data
cd OA
mkdir barcode_files count_matrices 
```

Script: ~/Downloads/RNAseq/memonet_github_repo/memonet/scripts/pairwiseOA_clean.r

Wd: ~/Downloads/RNAseq/AIBSmapping/OA

Query: MEMONET data
Reference: AIBS data

Input: 
- 10X directory of all mice combined: ~/Downloads/RNAseq/data/memonet_data/combined_cellranger_no-normalization/outs/filtered_feature_bc_matrix/
- AIBS dataset. The three files here: ~/Downloads/RNAseq/data/AIBS_data

Output: 
- prediction_scores.csv: lists each cell, the predicted cell type label, and prediction scores for each cell type
- umap_referenceAIBS.png: umap of the reference (AIBS data) with cell type labels
- umap_queryMEMONET.png: umap of the query (MEMONET data) projected onto AIBS space, labeled with the predicted labels

## Test the accuracy of mapping on AIBS data
### 1. Generate a test dataset (downsample to 25% of each cell type; remove sample cells from rest of reference) and perform label transfer from the remaining 75% of data. Do this 100 times.

Query: 25% of AIBS data

Reference: remaining 75% of AIBS data

Script: testA_prediction_cutoff.r

Wd: ~/Downloads/RNAseq/AIBSmapping/test

Input: 
- AIBS dataset. The three files here: RNAseq/data/AIBS_data
 
Output:
- 100 prediction_scores_*.csv files

### 2. Combine the 100 prediction score files into one file
Script: testA_prediction_cutoff_table.r

Wd: ~/Downloads/RNAseq/AIBSmapping/test

Input: 
- 100 prediction_scores_*.csv files
 
Output:
- prediction_cutoff.csv

### 3. Summarize results in RStudio: calculate false classification percentage, generate confusion matrix, compare mean and median scores of the test with OA mapping
Script: testA_prediction_cutoff2.r

Wd: ~/Downloads/RNAseq/AIBSmapping/test

Input: 
- prediction_cutoff.csv
- ~/Downloads/RNAseq/AIBSmapping/OA/prediction_scores.csv

Output:
- AIBStest_confusionMtx.png: heatmap of original cell type labels vs predicted labels for AIBS testing
- maxPredictionScores-AIBStest_and_OA.csv: table of mean and median prediction.score.max for AIBS testing and OA mapping; shows calculation for all celltypes and subset for L2/3 cells

## Choose L2/3 cutoff score of 0.3
Script: prediction_score_cutoff_barcodes.r
- This script generates a barcode list for various prediction score cutoff values. We chose a cutoff value of 0.3, meaning any cell that was predicted to be L2/3 and has a sum of prediction scores for L2/3 IT_1, L2/3 IT_2, L2/3 IT_3 <= 0.3 will not be included in downstream analysis.

Input: ~/Downloads/RNAseq/AIBSmapping/OA/prediction_scores.csv

Output: ~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/L23barcodes-fromAIBS_0.3.csv

# Investigate cell type proportions

Script: OverallDatasetDescription2.r, Part A

Functions:
- calculate cell subclass proportion (neurons and glia)
- calculate inhibitory neuron percentage
- calculate neuron subclass proportion (neurons only)
- plot neuron subclass proportion, comparing MEMONET dataset and AIBS

Input:
- AIBS metadata: ~/Downloads/RNAseq/data/AIBS_data/aibs_barcodes.tsv
- cell type predictions for MEMONET data: ~/Downloads/RNAseq/AIBSmapping/OA/prediction_scores.csv

## L2/3 subtypes: Is there train/control enrichment?
Script: statistics.r

Wd: ~/Downloads/RNAseq/AIBSmapping/OA/

Input: 
- Choose an input file at the beginning of script. For this part calculating L2/3 proportions, use ~/Downloads/RNAseq/AIBSmapping/OA/prediction_scores.csv
- L2/3 barcode list: ~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/L23barcodes-fromAIBS_0.3.csv

Output:
- L23subclass_tr_ctrl_prop_hist.png
- 'summary' variable lists p-values

# DE analysis for glutamatergic, GABAergic neurons
## 1. Generate barcode files for glutamatergic and GABAergic neurons
Script: OverallDatasetDescription2.r, Part B

Input: data_neurons variable from Part A

Output: 
- ~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/AIBS-defined_glut_barcodes.csv
- ~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/AIBS-defined_GABA_barcodes.csv

## 2. Run DESeq2, train vs control
```{r}
cd ~/Downloads/RNAseq/AIBSmapping/OA
mkdir DESeq2
cd DESeq2
mkdir glut_tr_vs_ctrl GABA_tr_vs_ctrl
```

Script: DESeq2_tr_vs_ctrl.r

**DE analysis of glutamatergic neurons:**

Wd: ~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/glut_tr_vs_ctrl

Input: ~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/AIBS-defined_glut_barcodes.csv

Output:
- 1unnormalized_counts_from_dds.csv: unnormalized gene expression
- 1normalized_sizeFactors_calculateSumFactors.csv: size factors that generate the normalized data
- 1normalized_counts_from_dds.csv: normalized gene expression
- 1_train_vs_control_all_genes.csv: DESeq2 results for all genes
- 1_train_vs_control_sig_genes.csv: DESeq2 results for significant genes (padj <0.05)

**DE analysis of GABAergic neurons:**

Wd: ~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/GABA_tr_vs_ctrl

Input: ~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/AIBS-defined_GABA_barcodes.csv

Output:
- 1unnormalized_counts_from_dds.csv: unnormalized gene expression
- 1normalized_sizeFactors_calculateSumFactors.csv: size factors that generate the normalized data
- 1normalized_counts_from_dds.csv: normalized gene expression
- 1_train_vs_control_all_genes.csv: DESeq2 results for all genes
- 1_train_vs_control_sig_genes.csv: DESeq2 results for significant genes (padj <0.05)

## 3. Summarize DE results: what IEGs are significant?
DEGvisuals.r (* update with heatmaps and lineplots)

Could also use genelist.r to make a table of which IEGs are up or down per cluster (* would have to format a version for manuscript)


# DE analysis for L2/3 neurons
```{r}
cd ~/Downloads/RNAseq/AIBSmapping/OA/DESeq2
mkdir L23_0.3_tr_vs_ctrl
```

## 1. Run DESeq2: train vs control
Script: DESeq2_tr_vs_ctrl.r

Wd: ~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/L23_0.3_tr_vs_ctrl

Input: L2/3 barcode list: ~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/L23barcodes-fromAIBS_0.3.csv

Output:
- 1unnormalized_counts_from_dds.csv: unnormalized gene expression
- 1normalized_sizeFactors_calculateSumFactors.csv: size factors that generate the normalized data
- 1normalized_counts_from_dds.csv: normalized gene expression
- 1_train_vs_control_all_genes.csv: DESeq2 results for all genes
- 1_train_vs_control_sig_genes.csv: DESeq2 results for significant genes (padj <0.05)

## 2. Summarize DE results: what IEGs are significant?
DEGvisuals.r

Could also use genelist.r to make a table of which IEGs are up or down per clusster (* would have to format a version for manuscript)

## 3. How many DEGs overlap with the 3000 experience-dependent genes (EDGs) used for clustering?
### a. Download EDG list
```{r}
cd ~/Downloads/RNAseq/
mkdir cluster_by_genes
```

Location on github:   *   ...../Updated_TopGenesAccordingtoLDA_trVsCtrl.csv

Download to: ~/Downloads/RNAseq/cluster_by_genes

### b. Calculate gene set overlap
Script: OverallDatasetDescription2.r, part C

Input:
- DG list: ~/Downloads/RNAseq/cluster_by_genes/Updated_TopGenesAccordingtoLDA_trVsCtrl.csv
- L2/3 DEG list: ~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/L23/all_cells_combined/1_train_vs_control_sig_genes.csv



# GO analysis of the 3000 EDGs
Before continuing with the analysis, we want to check whether the 3000 EDGs are actually related to learning.
## To run GO, you first need to generate a background gene list: this is all genes in the L2/3 dataset that are input to DESeq2
Script: background_genes-enrichment.r

Wd: ~/Downloads/RNAseq/AIBSmapping/OA/

Input: 
- L2/3 barcode list: ~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/L23barcodes-fromAIBS_0.2.csv

Output:
- background_genes-enrichment0.2.csv: background gene list, in gene symbol format
- background_genes-enrichment_entrezid0.2.csv: background gene list, in entrezid format

## Now run GO
```{r}
cd ~/Downloads/RNAseq/cluster_by_genes
mkdir 0.3cutoff
cd 0.3cutoff
mkdir GO
```

Script: Figures.r part “ GO analysis of classifier gene list”

Wd: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO

Input: 
- DG list: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/Updated_TopGenesAccordingtoLDA_trVsCtrl.csv
- Background gene list: ~/Downloads/RNAseq/AIBSmapping/OA/background_genes-enrichment0.2.csv

Output: 
- DiscriminantGenesGO_all.csv: GO results for the 1000 DGs (*update names)
- DiscriminantGenesGO_all.png: treeplot showing GO results for the 1000 DGs


Now that we've confirmed the EDGs do enrich for learning-related functions, we can cluster our L2/3 neurons using these genes.

# Cluster L2/3 neurons using EDGs
## 1. Normalize counts
Script: DESCnormalization.r

Wd: ~/Downloads/RNAseq/AIBSmapping/OA/count_matrices

Set 'cutoff' variable to '0.3'

Input:
- memonet data: ~/Downloads/RNAseq/data/memonet_data/combined_cellranger_no-normalization/outs/filtered_feature_bc_matrix

Output:
- ~/Downloads/RNAseq/QC/cells_after_QC.csv: list of cells remaining after QC steps on whole dataset. Needed for DESeq2 of clusters
- ~/Downloads/RNAseq/QC/genes_after_QC.csv: list of genes remaining after QC steps
- ~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/unnormalized_counts_L23_0.3.csv: unnormalized counts of L2/3 cells, 0.3 cutoff
- ~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/DESCnormalized_counts_L23_0.3.csv: normalized counts of L2/3 cells, 0.3 cutoff. This will be used as input for classifier training
- ~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/sampleIDs9.10.22.csv: list of cells and mouse ID, for reference during classifier training
- ~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/sampleIDs9.10.22_transpose.csv: same as above but the transposed version. Not sure which one Nathan used

## 2. Subset to EDGs 
Script: DESC_inputFormat.r

Wd: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff

Input:
- EDG list: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/PredictionGenesDescending.csv (*changed)
- L2/3 normalized counts: ~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/DESCnormalized_counts_L23_0.3.csv

Output: L23_0.3_DGmtx.csv

## 3. Run DESC clustering
```{r}
cd ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff
mkdir DESC
```

Script: DESC_DGclustering.py

Wd: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC

Input: 
- Normalized expression matrix, subset to DGs: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/L23_0.3_DGmtx.csv

Output:
- clusters.csv: lists each cell barcode and the cluster it is assigned
- umap.csv: umap x and y coordinates; row indices correspond to the barcode order in clusters.csv
- tsne.csv: tsne x and y coordinates
- result_DESC/: directory for encoder weights and model info
- figures/
   - umap0.8desc_0.8.png: visual of cluster umap projection 

## 4. Visualize cluster train/control proportion
Script: classifier_umap_plot.r

Wd: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/figures

Input: 
- Cluster file: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters.csv
- Umap coordinate file: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/umap.csv

Output: classifier_umap_tr_ctrl.svg

### Is there train/control enrichment?
Script: statistics.r

Wd: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/figures

Input:
- Choose an input file at the beginning of script. For this part calculating cluster train/control proportions, use ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters.csv

Output:
- cluster_tr_ctrl_prop_hist.svg
- 'summary' variable lists p-values

# DE analysis of clusters
## 1. Run DESeq2, one cluster vs the others
```{r}
cd ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC
mkdir DESeq2
```

Script: DESeq2_whole_subclusters.r

Wd: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2

Input: 
- Cluster file: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters.csv

Output directory: all_cells/
- unnormalized_counts_from_dds.csv: unnormalized gene expression
- normalized_sizeFactors_calculateSumFactors.csv: size factors that generate the normalized data
- normalized_counts_from_dds.csv: normalized gene expression
- *_vs_others_all_genes.csv: DESeq2 results for all genes
- *_vs_others_sig_genes.csv: DESeq2 results for significant genes (padj <0.05)

**Combine cluster result files into one file:** 

Script: DESeq2-table.r

Wd: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/all_cells

Input: files ending in *_all_genes.csv

Output: 
- DEGstats_allGenes.csv: all genes
- DEGstats_padj0.05.csv: significant genes (padj <0.05)

## 2. Summarize DE results: what IEGs are significant?
DEGvisuals.r

Could also use genelist.r to make a table of which IEGs are up or down per clusster (* would have to format a version for manuscript)

# Map AIBS L2/3 cells onto our L2/3 cells and annotate by cluster 
```{r}
cd ~/Downloads/RNAseq/AIBSmapping
mkdir AO
  # 'AO' stands for AIBS onto Our (memonet) cells
```

Script: labelTransfer_AO_umapDESC.r

Wd: ~/Downloads/RNAseq/AIBSmapping/AO

Input:
- Memonet data: ~/Downloads/RNAseq/data/memonet_data/combined_cellranger_no-normalization/outs/filtered_feature_bc_matrix/
- AIBS dataset. The three files here: ~/Downloads/RNAseq/data/AIBS_data
- Cluster file: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters.csv
- Umap coordinate file: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/umap.csv

Output:
- prediction_scores.csv: prediction score for the AIBS cells onto the clusters
- Images:
  - umap_classifierRef.png: our cells visualized in umap space 
  - umap_aibsL23_colored_by_classifierCluster.png: aibs cells placed in classifier space, colored by their predicted cluster

**What proportion of AIBS cells map to each cluster?**
Script: labelTransfer_AO_stats.r

Wd: ~/Downloads/RNAseq/AIBSmapping/AO

Input: *  ~/Downloads/RNAseq/AIBSmapping/AO/_/prediction_scores.csv

Output: *
- pie chart
- summary table

Since C_ receives the most AIBS cell assignments, we would like to assume C_ as the baseline condition and run another DE analysis of each cluster in reference to C_.

# DE analysis of clusters
## 1. Run DESeq2, one cluster vs baseline
Script: DESeq2_clusterX_vs_clusterY_whole.r

Wd: *  /work/pi_yingzhang_uri_edu/kdunton/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/

Input: 
- Cluster file: ~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters.csv

Output directory: * (the script makes a separate directory per cluster and then you combine the _all_genes files into a new dir for table script
- unnormalized_counts_from_dds.csv: unnormalized gene expression
- normalized_sizeFactors_calculateSumFactors.csv: size factors that generate the normalized data
- normalized_counts_from_dds.csv: normalized gene expression
- x_vs_4_all_genes.csv: DESeq2 results for all genes
- x_vs_4_sig_genes.csv: DESeq2 results for significant genes (padj <0.05)

**Combine cluster result files into one file:** 

Script: DESeq2-table.r

Wd: *

Input: files ending in *_all_genes.csv

Output: 
- DEGstats_allGenes.csv: all genes
- DEGstats_padj0.05.csv: significant genes (padj <0.05)

## 2. Summarize DE results: what IEGs are significant?
DEGvisuals.r

Could also use genelist.r to make a table of which IEGs are up or down per cluster (* would have to format a version for manuscript)







