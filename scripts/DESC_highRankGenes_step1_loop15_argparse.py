# this script is for running initial DESC clusters back through DESC to get subclusters
# also calculates silhouette score for cluster results


import desc as desc
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import episcanpy
import argparse

sc.settings.verbosity = 3  #verbosity: errors(0), warnings(1), info(2), hints(3)
sc.logging.print_versions()


#import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--neighbor',type=float)
parser.add_argument('--louvain',type=float)
args = parser.parse_args()


######## set parameters to loop through ########
#n_neighbors = [8,9,10,11,12,13,14,15,18,19,20,21,22,23,24,25]
#n_neighbors = [19,20,21,22,23,24,25]
#n_neighbors = np.array(n_neighbors)
n = args[0]
#louvain_resolution = [0.7,0.75,0.8,0.85,0.9]
#louvain_resolution = [0.4,0.45,0.5,0.55,0.6,0.65]
#louvain_resolution = np.array(louvain_resolution)
L = args[1]


######### run DESC #########
print("n_neighbors=" + str(n))
print("louvain_resolution=" + str(L))

path = '/work/pi_yingzhang_uri_edu/kdunton/RNAseq/cluster_by_genes/0.3cutoff/L23_0.3_DGmtx.csv'
print("Path:"+ str(path))
adata = sc.read(path)
print("adata initial")
print(adata)
adata.raw = adata

## z-score the genes
desc.scale(adata, zero_center=True, max_value=3)

save_dir = ["result_DESC22.n"+str(n)+".L"+str(L)]
#save_dir is a list
#convert to string
save_dir = ''.join(map(str,save_dir))

adata = desc.train(adata, dims=[adata.shape[1], 64, 32], tol=0.001, n_neighbors=n,
                   batch_size=256, louvain_resolution=L,
                   save_dir=save_dir, do_tsne=True, learning_rate=300,
                   do_umap=True, num_Cores_tsne=4,
                   save_encoder_weights=False)


######### Visualize #########
#Umap: plot the clusterIDs
basis = ["umap" + str(L)]
basis = ''.join(map(str,basis))
color = ["desc_" + str(L)]
color = ''.join(map(str,color))
save = ["desc.dec22.n" + str(n) + ".L" + str(L) + ".png"]
print(save)
save = ''.join(map(str,save))

sc.pl.scatter(adata, basis=basis, color=color, save = save)



######### save result #########
# umap coordinates:
umap = ["X_umap" + str(L)]
umap = ''.join(map(str,umap))
obsm_data=pd.DataFrame(adata.obsm[umap])
save_umap = ["umap_dec22.n" + str(n) + ".L" + str(L) + ".csv"]
save_umap = ''.join(map(str,save_umap))
obsm_data.to_csv(save_umap,sep=",")

# clustering result/labels:
clu = ["desc_" + str(L)]
clu = ''.join(map(str,clu))
obs_data=pd.DataFrame(adata.obs[clu])
save_clu = ["clusters_dec22.n" + str(n) + ".L" + str(L) + ".csv"]
save_clu = ''.join(map(str,save_clu))
obs_data.to_csv(save_clu,sep=",")



######### calculate silhouette score #########
cluster_annot = ['desc_' + str(L)]
cluster_annot = ''.join(map(str,cluster_annot))
value = ['X_Embeded_z' + str(L)]
value = ''.join(map(str,value))
episcanpy.tl.silhouette(adata, cluster_annot=cluster_annot, value=value)
# saves individual silhouette scores to obs, with cluster IDs

sil = pd.DataFrame(adata.obs)
save_sil = ["sil_scores_n" + str(n) + ".L" + str(L) + ".csv"]
save_sil = ''.join(map(str,save_sil))
sil.to_csv(save_sil,sep=',')
print("************ saved sil ***************")
print(save_sil)


print("done!")