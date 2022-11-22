## this script summarizes the number of positive and negative reactivation related genes in MEMONET clusters
# calculates reactivation score
# makes plots of reactivation score counts 


library(tidyverse)
library(schoolmath) #for is.positive
library(RColorBrewer)



### read in DEG data
# each cluster vs others as ref
genes <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/all_cells/DEGstats_padj0.05.csv")

#-- reactivation related genes, ie reference list  
genelist_reactivated <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/all_cells/Jaeger_DG_React_NotReact.csv")
genelist_reactivated$boolean <- is.positive(genelist_reactivated$test_statistic)
genelist_reactivated$boolean <- as.character(genelist_reactivated$boolean)
  #need to convert from logical to character type for recode to work
genelist_reactivated <- genelist_reactivated %>% mutate(result = recode(boolean,
                                                                          'TRUE'='upregulated',
                                                                          'FALSE'='downregulated'))
genelist_reactivated_up <- subset(genelist_reactivated, result %in% 'upregulated')
genelist_reactivated_up <- genelist_reactivated_up[,'genesymbol']
genelist_reactivated_down <- subset(genelist_reactivated, result %in% 'downregulated')
genelist_reactivated_down <- genelist_reactivated_down[,'genesymbol']
#--

# subset genes to the reference list
genes <- subset(genes, padj < 0.05)
genes_upregReact <- subset(genes, gene %in% genelist_reactivated_up)
genes_downregReact <- subset(genes, gene %in% genelist_reactivated_down)


print("Positive reactivation genes: how many are up- and down-regulated in our clusters?")
table(genes_upregReact$cluster, genes_upregReact$result)
table_upregReact <- as.data.frame(table(genes_upregReact$cluster, genes_upregReact$result))
print("Negative reactivation genes: how many are up- and down-regulated in our clusters?")
table(genes_downregReact$cluster, genes_downregReact$result)
table_downregReact <- as.data.frame(table(genes_downregReact$cluster, genes_downregReact$result))


## calculate the strength (we'll call it 'reactivation score') of consistent directional regulation among the reactivation-related genes
clusters <- sort(unique(genes$cluster))
summary_reactivation <- data.frame(matrix(ncol = 6, nrow = length(clusters)))
colnames(summary_reactivation) <- c('cluster','reactivation_score','up.up','down.down','up.down','down.up')
summary_reactivation$cluster <- clusters

for(cluster in clusters){
  # first deal with the upregulated reactivation genes
  x_up <- subset(table_upregReact, Var1 %in% cluster)
  up.up <- x_up$Freq[x_up$Var2 == 'upregulated' ]
  #this is the number of genes upregulated in the cluster and also upregulated in the reference
  down.up <- x_up$Freq[x_up$Var2 == 'downregulated' ]
  #this is the number of genes downregulated in the cluster but upregulated in the reference 
  
  # now the downregulated reactivation genes
  x_down <- subset(table_downregReact, Var1 %in% cluster)
  down.down <- x_down$Freq[x_down$Var2 == 'downregulated' ]
  #this is the number of genes downregulated in the cluster and also downregulated in the reference 
  up.down <- x_down$Freq[x_down$Var2 == 'upregulated' ]
  #this is the number of genes upregulated in the cluster but downregulated in the reference 
  
  ## calculate reactivation score
  reactivation_score <- ((up.up + down.down) / (down.up + up.down)) - 1 # subtract 1 at the end so that the chance level is at 0 instead of 1
  ## append to df
  summary_reactivation$reactivation_score[summary_reactivation$cluster == cluster] <- reactivation_score
  summary_reactivation$up.up[summary_reactivation$cluster == cluster] <- up.up
  summary_reactivation$down.down[summary_reactivation$cluster == cluster] <- down.down
  summary_reactivation$up.down[summary_reactivation$cluster == cluster] <- up.down
  summary_reactivation$down.up[summary_reactivation$cluster == cluster] <- down.up
}
write.csv(summary_reactivation, "~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/all_cells/ReactivationScores.csv", row.names=FALSE)
# nomenclature of 'down.up' variable, for example: 
  # before the period means the sign for the cluster, ie downregulated in the cluster. 
  # After the period means the sign for the reference reactivated gene list, ie upregulated in the list
  # so this is the number of genes that are downregulated in the cluster but upregulated in the reference list



#### plot
summary_reactivation$cluster <- as.factor(summary_reactivation$cluster)

## reactivation score:
#line plot:
ggplot(summary_reactivation, aes(x=cluster, y=reactivation_score, group=1)) +
  #geom_line() +
  geom_point(fill='black',shape=23) +
  scale_y_continuous(limits=c(-0.5,2)) +
  geom_hline(yintercept=0, linetype=2, col = 'black') +
  theme_classic() +
  xlab('Cluster')+ ylab('Reactivation score')
ggsave("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/all_cells/reactivation_line.svg")


## plot the number of genes for pro-reactivation and anti-reactivation per cluster
# make the pro-reactivation gene (up.up/down.down) counts positive; make anti-reactivation gene (up.down/down.up) counts negative
# make a 'direction' column that's TRUE for up.up/down.down and FALSE for up.down/down.up 
  #https://stackoverflow.com/questions/52906237/how-to-make-column-values-negative-if-false-in-another-column

# pivot table
summary_reactivation2 <- summary_reactivation
summary_reactivation2$reactivation_score <- NULL
summary_reactivation2_long <- gather(summary_reactivation2, key=match, value=count, up.up:down.up, factor_key=TRUE)
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)

summary_reactivation2_long <- summary_reactivation2_long %>% 
  mutate(direction = case_when(match == 'up.up' ~ "TRUE", match == 'down.down' ~ "TRUE",
                               match == 'up.down' ~ "FALSE", match == 'down.up' ~ "FALSE"))
  # when 'direction' is TRUE, the 'count' column value stays positive; when FALSE, the count column value becomes negative 
summary_reactivation2_long <- summary_reactivation2_long %>%
  mutate(count_direction = ifelse(direction, count, -count))

ggplot(summary_reactivation2_long, aes(x=cluster, y=count_direction, fill=match)) +
  geom_bar(stat='identity', position='stack', width=0.5) +
  scale_fill_manual(values=c('red','darksalmon','blue','lightblue')) +
  geom_hline(yintercept=0, linetype=1, col = 'black') +
  theme_classic() +
  xlab('Cluster')+ ylab('Number of genes')
ggsave("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/all_cells/reactivation_matches.svg")

 
