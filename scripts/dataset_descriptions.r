# this script calculates various things:
# Part A: cell type composition comparison of AIBS and MEMONET datasets
#         inhibitory neuron percentage of AIBS and MEMONET datasets
# Part B: save barcode files for specified neuron types
#         proportion of AIBS and MEMONET cells in each L2/3 subclass (not split by train/control)
# Part C: gene overlap between L2/3 DE analysis and the 3000 EDGs
# Part D: percentage of MEMONET cells in each L2/3 cluster
# Part E: proportion of L2/3 AIBS cells mapped to L2/3 clusters


library(tidyverse)
library(RColorBrewer) #for colorblind palette 
library(matrixStats)
library(svglite)


#-------- Part A: Cell type composition of AIBS and MEMONET datasets; inhibitory percentages --------# 
#### first get aibs proportions ####
aibs <- read.csv("~/Downloads/RNAseq/data/AIBS_data/aibs_barcodes.tsv") 
aibs <- aibs[,c('X','cluster_label')]
aibs$AIBS_predicted_label <- aibs$cluster_label
aibs$cluster_label <- NULL


## glut GABA glial subclass percentages ##
# condense the specific labels together into subclasses
aibs_data <- aibs %>%
  mutate(subclass = recode(AIBS_predicted_label,
                           "L2/3 IT_1" = "L2/3 IT",
                           "L2/3 IT_2" = "L2/3 IT",
                           "L2/3 IT_3" = "L2/3 IT",
                           "L5 ET_1" = "L5 ET",
                           "L5 ET_2" = "L5 ET",
                           "L5 ET_3" = "L5 ET",
                           "L5 ET_4" = "L5 ET",
                           "L5 ET_5" = "L5 ET",
                           "L5 ET_6" = "L5 ET",
                           "L5 IT Pld5_1" = "L5 IT",
                           "L5 IT Pld5_2" = "L5 IT",
                           "L5 IT Rspo1_1" = "L5 IT",
                           "L5 IT Rspo1_2" = "L5 IT",
                           "L5 IT Rspo1_3" = "L5 IT",
                           "L5 IT S100b" = "L5 IT",
                           "L5 NP Slc17a8_1" = "L5 NP",
                           "L5 NP Slc17a8_2" = "L5 NP",
                           "L6 CT Brinp3" = "L6 CT",
                           "L6 CT Cpa6" = "L6 CT",
                           "L6 CT Gpr139" = "L6 CT",
                           "L6 CT Nxph2" = "L6 CT",
                           "L6 IT Sulf1_1" = "L6 IT",
                           "L6 IT Sulf1_2" = "L6 IT",
                           "L6 IT Sulf1_3" = "L6 IT",
                           "L6 IT Sulf1_4" = "L6 IT",
                           "L6 NP CT" = "L6 NP",
                           "L6 NP Trh_1" = "L6 NP",
                           "L6 NP Trh_2" = "L6 NP",
                           "L6b Kcnip1" = "L6b",
                           "L6b Ror1" = "L6b",
                           "L6b Rprm_1" = "L6b",
                           "L6b Rprm_2" = "L6b",
                           "L6b Shisa6" = "L6b",
                           "Lamp5 Egln3" = "Lamp5",
                           "Lamp5 Lhx6" = "Lamp5",
                           "Lamp5 Pax6" = "Lamp5",
                           "Lamp5 Pdlim5" = "Lamp5",
                           "Lamp5 Slc35d3" = "Lamp5",
                           "Pvalb Gabrg1" = "Pvalb",
                           "Pvalb Gpr149" = "Pvalb",
                           "Pvalb Il1rapl2_1" = "Pvalb",
                           "Pvalb Il1rapl2_2" = "Pvalb",
                           "Pvalb Prdm8" = "Pvalb",
                           "Pvalb Reln" = "Pvalb",
                           "Pvalb Vipr2" = "Pvalb",
                           "Sncg Col14a1" = "Sncg",
                           "Sncg Slc17a8" = "Sncg",
                           "Sst C1ql3_1" = "Sst",
                           "Sst C1ql3_2" = "Sst",
                           "Sst Calb2" = "Sst",
                           "Sst Chodl" = "Sst",
                           "Sst Hpse" = "Sst",
                           "Sst Myh8" = "Sst",
                           "Sst Pvalb Etv1_1" = "Sst",
                           "Sst Pvalb Etv1_2" = "Sst",
                           "Sst Pvalb Etv1_3" = "Sst",
                           "Vip Chat" = "Vip",
                           "Vip Crispld2" = "Vip",
                           "Vip Htr1f" = "Vip",
                           "Vip Igfbp6" = "Vip",
                           "Vip Serpinf1" = "Vip",
                           "Astro Aqp4" = "Astrocyte",
                           "Endo Slc38a5" = "Endothelial",
                           "Macrophage" = "Microglia",
                           "Oligo Opalin" = "Oligodendrocyte",
                           "OPC Pdgfra_1" = "OPC",
                           "OPC Pdgfra_2" = "OPC"))
table_subclass_aibs <- as.data.frame(table(aibs_data$subclass))
table_subclass_aibs$sum <- sum(table_subclass_aibs$Freq)
table_subclass_aibs$prop <- table_subclass_aibs$Freq / table_subclass_aibs$sum
# add dataset column so this df can be compared with the one generated in next section for MEMONET dataset, in order to plot both datasets
table_subclass_aibs$dataset <- 'AIBS'


## calculate inhibitory neuron percentage: GABA count / total neuron count
aibs_data <- aibs_data %>%
  mutate(class = recode(subclass,
                        "L2/3 IT" = "glut",
                        "L5 ET" = "glut",
                        "L5 IT" = "glut",
                        "L5 NP" = "glut",
                        "L6 CT" = "glut",
                        "L6 IT" = "glut",
                        "L6 NP" = "glut",
                        "L6b" = "glut",
                        "Lamp5" = "GABA",
                        "Pvalb" = "GABA",
                        "Sncg" = "GABA",
                        "Sst" = "GABA",
                        "Vip" = "GABA",
                        "Astrocyte" = "glial",
                        "Endothelial" = "glial",
                        "Microglia" = "glial",
                        "Oligodendrocyte" = "glial",
                        "OPC" = "glial"))
table_broad_aibs <- as.data.frame(table(aibs_data$class))

neuron_total_aibs <- sum(c(table_broad_aibs$Freq[table_broad_aibs$Var1 == 'GABA'], table_broad_aibs$Freq[table_broad_aibs$Var1 == 'glut']))
GABA_total_aibs <- table_broad_aibs$Freq[table_broad_aibs$Var1 == 'GABA']
inhib_percent_aibs <- (GABA_total_aibs / neuron_total_aibs)*100
print("inhibitory neuron percentage, AIBS dataset:")
inhib_percent_aibs


## subclass proportion of just neurons 
aibs_data_neurons <- subset(aibs_data, class %in% c('glut','GABA'))
table_neurons_aibs <- as.data.frame(table(aibs_data_neurons$subclass))
table_neurons_aibs$sum <- sum(table_neurons_aibs$Freq)
table_neurons_aibs$prop <- table_neurons_aibs$Freq / table_neurons_aibs$sum
table_neurons_aibs$dataset <- 'AIBS'





#### now get MEMONET dataset proportions, from OA mapping ####
datai <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/prediction_scores.csv")
datai <- datai[,c('X','AIBS_predicted_label')]


## glut GABA glial subclass percentages ##
# condense the specific labels together into subclasses
data <- datai %>%
  mutate(subclass = recode(AIBS_predicted_label,
                            "L2/3 IT_1" = "L2/3 IT",
                            "L2/3 IT_2" = "L2/3 IT",
                            "L2/3 IT_3" = "L2/3 IT",
                            "L5 ET_1" = "L5 ET",
                            "L5 ET_2" = "L5 ET",
                            "L5 ET_3" = "L5 ET",
                            "L5 ET_4" = "L5 ET",
                            "L5 ET_5" = "L5 ET",
                            "L5 ET_6" = "L5 ET",
                            "L5 IT Pld5_1" = "L5 IT",
                            "L5 IT Pld5_2" = "L5 IT",
                            "L5 IT Rspo1_1" = "L5 IT",
                            "L5 IT Rspo1_2" = "L5 IT",
                            "L5 IT Rspo1_3" = "L5 IT",
                            "L5 IT S100b" = "L5 IT",
                            "L5 NP Slc17a8_1" = "L5 NP",
                            "L5 NP Slc17a8_2" = "L5 NP",
                            "L6 CT Brinp3" = "L6 CT",
                            "L6 CT Cpa6" = "L6 CT",
                            "L6 CT Gpr139" = "L6 CT",
                            "L6 CT Nxph2" = "L6 CT",
                            "L6 IT Sulf1_1" = "L6 IT",
                            "L6 IT Sulf1_2" = "L6 IT",
                            "L6 IT Sulf1_3" = "L6 IT",
                            "L6 IT Sulf1_4" = "L6 IT",
                            "L6 NP CT" = "L6 NP",
                            "L6 NP Trh_1" = "L6 NP",
                            "L6 NP Trh_2" = "L6 NP",
                            "L6b Kcnip1" = "L6b",
                            "L6b Ror1" = "L6b",
                            "L6b Rprm_1" = "L6b",
                            "L6b Rprm_2" = "L6b",
                            "L6b Shisa6" = "L6b",
                            "Lamp5 Egln3" = "Lamp5",
                            "Lamp5 Lhx6" = "Lamp5",
                            "Lamp5 Pax6" = "Lamp5",
                            "Lamp5 Pdlim5" = "Lamp5",
                            "Lamp5 Slc35d3" = "Lamp5",
                            "Pvalb Gabrg1" = "Pvalb",
                            "Pvalb Gpr149" = "Pvalb",
                            "Pvalb Il1rapl2_1" = "Pvalb",
                            "Pvalb Il1rapl2_2" = "Pvalb",
                            "Pvalb Prdm8" = "Pvalb",
                            "Pvalb Reln" = "Pvalb",
                            "Pvalb Vipr2" = "Pvalb",
                            "Sncg Col14a1" = "Sncg",
                            "Sncg Slc17a8" = "Sncg",
                            "Sst C1ql3_1" = "Sst",
                            "Sst C1ql3_2" = "Sst",
                            "Sst Calb2" = "Sst",
                            "Sst Chodl" = "Sst",
                            "Sst Hpse" = "Sst",
                            "Sst Myh8" = "Sst",
                            "Sst Pvalb Etv1_1" = "Sst",
                            "Sst Pvalb Etv1_2" = "Sst",
                            "Sst Pvalb Etv1_3" = "Sst",
                            "Vip Chat" = "Vip",
                            "Vip Crispld2" = "Vip",
                            "Vip Htr1f" = "Vip",
                            "Vip Igfbp6" = "Vip",
                            "Vip Serpinf1" = "Vip",
                            "Astro Aqp4" = "Astrocyte",
                            "Endo Slc38a5" = "Endothelial",
                            "Macrophage" = "Microglia",
                            "Oligo Opalin" = "Oligodendrocyte",
                            "OPC Pdgfra_1" = "OPC",
                            "OPC Pdgfra_2" = "OPC"))
table_subclass <- as.data.frame(table(data$subclass))
table_subclass$sum <- sum(table_subclass$Freq)
table_subclass$prop <- table_subclass$Freq / table_subclass$sum
table_subclass$dataset <- 'MEMONET'


## calculate inhibitory neuron percentage: GABA count / total neuron count
data <- data %>%
  mutate(class = recode(subclass,
                            "L2/3 IT" = "glut",
                            "L5 ET" = "glut",
                            "L5 IT" = "glut",
                            "L5 NP" = "glut",
                            "L6 CT" = "glut",
                            "L6 IT" = "glut",
                            "L6 NP" = "glut",
                            "L6b" = "glut",
                            "Lamp5" = "GABA",
                           "Pvalb" = "GABA",
                           "Sncg" = "GABA",
                           "Sst" = "GABA",
                           "Vip" = "GABA",
                           "Astrocyte" = "glial",
                           "Endothelial" = "glial",
                           "Microglia" = "glial",
                           "Oligodendrocyte" = "glial",
                           "OPC" = "glial"))
table_broad <- as.data.frame(table(data$class))

neuron_total <- sum(c(table_broad$Freq[table_broad$Var1 == 'GABA'], table_broad$Freq[table_broad$Var1 == 'glut']))
GABA_total <- table_broad$Freq[table_broad$Var1 == 'GABA']
inhib_percent <- (GABA_total / neuron_total)*100
print("inhibitory neuron percentage, MEMONET dataset:")
inhib_percent


## subclass proportion of just neurons 
data_neurons <- subset(data, class %in% c('glut','GABA'))
table_neurons <- as.data.frame(table(data_neurons$subclass))
table_neurons$sum <- sum(table_neurons$Freq)
table_neurons$prop <- table_neurons$Freq / table_neurons$sum
table_neurons$dataset <- 'MEMONET'




### plot neuron subclass proportion by dataset
table_bothDatasets_neurons <- rbind(table_neurons, table_neurons_aibs)
table_bothDatasets_neurons <- table_bothDatasets_neurons %>% rename(subclass=Var1)
# check that the sum of proportions per dataset is 1
check_table <- table_bothDatasets_neurons %>%
  group_by(dataset) %>%
  summarise(Sum = sum(prop))


# plot
colors_colorblind <- brewer.pal(n = 12, name = "Paired")
colors <- c(colors_colorblind,'gray')
ggplot(data=table_bothDatasets_neurons, aes(x=dataset, y=prop, fill=subclass)) +
  geom_bar(stat='identity', position='stack', width=0.8) +
  scale_fill_manual(values=colors) +
  xlab('Neuronal cell type') + ylab('Proportion per dataset') +
  theme_classic()
ggsave("~/Downloads/RNAseq/AIBSmapping/OA/CelltypePropPerDataset.svg")
ggsave("~/Downloads/RNAseq/AIBSmapping/OA/CelltypePropPerDataset.png")






#-------- Part B: Save file of neuron subtype barcodes for DE analysis --------#
### save file of glut barcodes, for DE analysis 
data_glut <- subset(data_neurons, class %in% 'glut')
data_glut <- as.data.frame(data_glut[,'X'])
colnames(data_glut) <- 'x'
write.csv(data_glut, "~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/AIBS-defined_glut_barcodes.csv", row.names=FALSE)

### save file of GABA barcodes, for DE analysis
data_GABA <- subset(data_neurons, class %in% 'GABA')
data_GABA <- as.data.frame(data_GABA[,'X'])
colnames(data_GABA) <- 'x'
write.csv(data_GABA, "~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/AIBS-defined_GABA_barcodes.csv", row.names=FALSE)

### save file of L5 neuron barcodes for reviewer analysis
data_L5 <- subset(data_neurons, subclass %in% c('L5 IT', 'L5 ET', 'L5 NP'))
data_L5 <- as.data.frame(data_L5[,'X'])
colnames(data_L5) <- 'x'
write.csv(data_L5, "~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/L5_barcodes.csv", row.names=FALSE)

## save file of GABA vip
data_vip <- subset(data_neurons, subclass %in% 'Vip')
data_vip <- as.data.frame(data_vip[,'X'])
colnames(data_vip) <- 'x'
write.csv(data_vip, "~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/vip_barcodes.csv", row.names=FALSE)

## save file of GABA sst
data_sst <- subset(data_neurons, subclass %in% 'Sst')
data_sst <- as.data.frame(data_sst[,'X'])
colnames(data_sst) <- 'x'
write.csv(data_sst, "~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/sst_barcodes.csv", row.names=FALSE)

## save file of GABA pvalb
data_pvalb <- subset(data_neurons, subclass %in% 'Pvalb')
data_pvalb <- as.data.frame(data_pvalb[,'X'])
colnames(data_pvalb) <- 'x'
write.csv(data_pvalb, "~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/pvalb_barcodes.csv", row.names=FALSE)





#-------- Part C: how many DEG from all L2/3 train vs control overlap with the 3000 EDGs? --------#
# read in classifier gene list
genelist <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/PredictionGenesDescending0.3.csv", header=FALSE)
genelist <- na.omit(genelist)
# subset to top 3000 (EDGs)
genelist <- genelist[1:3000,]

# read in DEG list
DEG <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/L23_0.3_tr_vs_ctrl/train_vs_control_sig_genes.csv")

#see what genes are the same from DEG list to EDG list
DEGoverlap <- subset(DEG, gene %in% genelist)
print("number of shared genes:")
length(DEGoverlap$gene)



#-------- Part D: percentage of cells in each cluster--------#
clusters <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters_n25.L0.65.csv")
names(clusters)[2] <- 'cluster'
table <- as.data.frame(table(clusters$cluster))
table$percent <- (table$Freq / length(clusters$X) )*100
# check the percent column adds to 100
sum(table$percent)

# visualize in pie chart
colors <- brewer.pal(n = 6, name = "Set2")
# to display the values, you have to write them in the label
pie_labels <- paste0(table$Var1, " = ", round(table$percent,2), "%")
pie(x=table$Freq, labels=pie_labels, clockwise=TRUE, col=colors, #labels=table$cluster
    main='Proportion of MEMONET cells per cluster')

names(table)[1] <- 'cluster'
write.csv(table, "~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/figures/cluster_percentage_tbl.csv", row.names = FALSE)





#-------- Part E: AO mapping --------#
# visualize proportion of AIBS cells mapped to each L2/3 cluster
#### read in prediction data ####
data <- read.csv("~/Downloads/RNAseq/AIBSmapping/AO/parameter_loop/n25.L0.65_prediction_scores.csv")
data <- data[,c('cluster_label','predicted.id')]


#### get observed counts of AIBS cells per classifier cluster ####
table <- as.data.frame(table(data$predicted.id))
colnames(table) <- c('cluster','Freq')
table$percent <- 100 * table$Freq/sum(table$Freq)
write.csv(table, "~/Downloads/RNAseq/AIBSmapping/AO/n25.L0.65/AIBSmap_table.csv", row.names = FALSE)
table$percent <- round(100 * table$Freq/sum(table$Freq), 2)

## pie chart
#https://r-coder.com/pie-chart-r/#Pie_chart_in_R_with_percentage
# set colors using a color palette 
colors <- brewer.pal(n = 8, name = "Set2")
colors <- c(colors, "gray")
# to display the values, you have to write them in the label
pie_labels <- paste0(table$cluster, " : ", table$percent, "%")
svg(filename="~/Downloads/RNAseq/AIBSmapping/AO/n25.L0.65/AIBSpie.svg")
pie(x=table$Freq, labels=pie_labels, clockwise=TRUE, col=colors,
    main='Proportion of AIBS cells')
dev.off()





