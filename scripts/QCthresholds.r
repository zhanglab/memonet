# Visualize data to determine proper QC thresholds


library(dplyr)
library(Seurat)
#library(patchwork) may not need this

options(future.globals.maxSize = 7000 * 1024^2)


#############  Load the data of each mouse ############# 
M262.data = Read10X(data.dir = "~/Downloads/RNAseq/data/memonet_data/data_download/slPsiwmg_JB_262_1_2_3/filtered_feature_bc_matrix/")
M262 = CreateSeuratObject(counts = M262.data, min.cells = 3, project = "M262")

M263.data = Read10X(data.dir = "~/Downloads/RNAseq/data/memonet_data/data_download/QF9GKKgY_JB_263_1_2_3/filtered_feature_bc_matrix/")
M263 = CreateSeuratObject(counts = M263.data, min.cells = 3, project = "M263")

M276.data = Read10X(data.dir = "~/Downloads/RNAseq/data/memonet_data/data_download/JB_276_1_2_3/outs/filtered_feature_bc_matrix/")
M276 = CreateSeuratObject(counts = M276.data, min.cells = 3, project = "M276")

M277.data = Read10X(data.dir = "~/Downloads/RNAseq/data/memonet_data/data_download/JB_277_1_2_3/outs/filtered_feature_bc_matrix/")
M277 = CreateSeuratObject(counts = M277.data, min.cells = 3, project = "M277")

M278.data = Read10X(data.dir = "~/Downloads/RNAseq/data/memonet_data/data_download/JB_278_1_2_3/outs/filtered_feature_bc_matrix/")
M278 = CreateSeuratObject(counts = M278.data, min.cells = 3, project = "M278")

M279.data = Read10X(data.dir = "~/Downloads/RNAseq/data/memonet_data/data_download/FHuOVucK_JB_279_1_2/filtered_feature_bc_matrix/")
M279 = CreateSeuratObject(counts = M279.data, min.cells = 3, project = "M279")


##############  Merging control and trained animals ############# 
ctrl = merge(M262, y = c(M276,M278), add.cell.ids = c("C1", "C2", "C3"), project = "C_mice")
ctrl$stim <- "CTRL"

train = merge(M263, y = c(M277,M279), add.cell.ids = c("T1", "T2", "T3"), project = "T_mice")
train$stim <- "TRAIN"


##### calculate mitochondiral gene ratio #####
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^mt-")
train[["percent.mt"]] <- PercentageFeatureSet(train, pattern = "^mt-")


# Visualize QC metrics as a violin plot
bitmap("ctrl_without_cutoff.png", width = 13, height = 11, units = 'in', res = 300)
VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

bitmap("train_without_cutoff.png", width = 13, height = 11, units = 'in', res = 300)
VlnPlot(train, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()


######## visualize after QC- once you look at above plots, decide on thresholds here and view how it cuts the data ########
ctrl <- subset(ctrl, subset = nFeature_RNA > 200 & nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt <1)
train <- subset(train, subset = nFeature_RNA > 200 & nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt <1)

bitmap("QC_ctrl.png", width = 13, height = 11, units = 'in', res = 300)
VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()
bitmap("QC_train.png", width = 13, height = 11, units = 'in', res = 300)
VlnPlot(train, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

