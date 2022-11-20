# this script is for running DESC using only a subset of genes (the 3000 EDGs)


import desc as desc
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import episcanpy
sc.settings.verbosity = 3  #verbosity: errors(0), warnings(1), info(2), hints(3)
sc.logging.print_versions()



######### load andata #########
path = '~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/L23_0.3_EDGmtx.csv'
print("Path:"+ str(path))
adata = sc.read(path)
  #this reads a counts matrix into anndata object; requires cells in rownames and genes in columns


######### preprocessing #########
#the counts in adata have already been cell-level and log normalized, so skip normalize_by_cell() and logp1()

## view object
print("adata initial")
print(adata)

#Set the .raw attribute of AnnData object to the logarithmized raw gene expression for downstream analysis,
# such as differential expression analysis, and pseudotime analysis.
# This simply freezes the state of the current AnnData object
adata.raw = adata

## z-score the genes
desc.scale(adata, zero_center=True, max_value=3)



######### run DESC #########
## set parameters
n = 25
print("n_neighbors:")
print(n)
res = 0.65
print("louvain_resolution:")
print(res)

save_dir = ["result_DESC.n"+str(n)+".L"+str(res)]
#save_dir is a list, convert to string
save_dir = ''.join(map(str,save_dir))

adata = desc.train(adata, dims=[adata.shape[1], 64, 32], tol=0.001, n_neighbors=n,
                   batch_size=256, louvain_resolution=[res],
                   save_dir=save_dir, do_tsne=True, learning_rate=300,
                   do_umap=True, num_Cores_tsne=4,
                   save_encoder_weights=True)




######### Visualize #########
#Umap plot the clusterIDs
basis = ["umap" + str(res)]
basis = ''.join(map(str,basis))
print(basis)
color = ["desc_" + str(res)]
color = ''.join(map(str,color))
print(color)
save = ["desc.n" + str(n) + ".L" + str(res) + ".png"]
print(save)
save = ''.join(map(str,save))
sc.pl.scatter(adata, basis=basis, color=color, save = save)


######### save result #########
## save AnnData object into .h5ad file ##
save_adata = ["desc_result.n" + str(n) + ".L" + str(res) + ".h5ad"]
save_adata = ''.join(map(str,save_adata))
adata.write(save_adata)

## save into csv ##
# umap coordinates:
umap = ["X_umap" + str(res)]
umap = ''.join(map(str,umap))
obsm_data=pd.DataFrame(adata.obsm[umap])
save_umap = ["umap_n" + str(n) + ".L" + str(res) + ".csv"]
save_umap = ''.join(map(str,save_umap))
obsm_data.to_csv(save_umap,sep=",")
# clustering result/labels:
clu = ["desc_" + str(res)]
clu = ''.join(map(str,clu))
obs_data=pd.DataFrame(adata.obs[clu])
save_clu = ["clusters_n" + str(n) + ".L" + str(res) + ".csv"]
save_clu = ''.join(map(str,save_clu))
obs_data.to_csv(save_clu,sep=",")


print("done!")
