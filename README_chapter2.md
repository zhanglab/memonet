This file details all steps of the analysis for Chapter 2, including scripts and location of input/output data files.



# Cluster by cell type with DESC
1. Determine QC thresholds through metric visualization
- Script: /data/zhanglab/kdunton/neuron_model/katie-scripts/RNAseq/scripts_unorganized/test_QCthresholdsDESC.r 
- Usage: sbatch /data/zhanglab/kdunton/neuron_model/katie-scripts/RNAseq/scripts_unorganized/test_QCthresholdsDESC.sh
- Wd: /data/zhanglab/kdunton/6samples_cluster/deepseq_3_clustering/DESC/DESC_final_run/figures
- Input: location of each mouse folder, ie /data/zhanglab/jingwang/brain/RNAseq/Takaki/Deep/deepseq_2/slPsiwmg_JB_262_1_2_3/filtered_feature_bc_matrix/
- Output: 
  - ctrl_without_cutoff.png
  - QC_ctrl.png (or ctrl_with_cutoff30000.png)
  - train_without_cutoff.png
  - QC_train.png (or train_with_cutoff30000.png)

2. Cluster the cells
- Script: /data/zhanglab/kdunton/neuron_model/katie-scripts/RNAseq/DESC_step1.py
- Usage: sbatch /data/zhanglab/kdunton/neuron_model/katie-scripts/RNAseq/DESC_step1.sh
- Wd: /data/zhanglab/kdunton/6samples_cluster/deepseq_3_clustering/DESC/DESC_final_run
- Input:
  - Directory to 10X cellranger aggr folder: /data/zhanglab/kdunton/6samples_cluster/deepseq_3_clustering/scrattch/combined_cellranger_no-normalization/outs/filtered_feature_bc_matrix
  - Inside filtered_feature_bc_matrix you need to unzip the 3 files (matrix, barcodes, features) and rename features to cells.tsv
- Output: 
  - Slurm output file will give the number of cells per cluster (near bottom)
  - clusters.csv
  - umap.csv
  - tsne.csv
  - figures/umap0.8desc_0.8.png - umap with each cluster a different color
  - result_DESC/ - contains files on the encoders and models
  - desc_allgenes.h5ad - anndata object containing all genes passing QC
  - desc_result.h5ad- anndata object containing only HVG genes

