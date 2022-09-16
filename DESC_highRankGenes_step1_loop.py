# this script is for running initial DESC clusters back through DESC to get subclusters


import desc as desc
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
sc.settings.verbosity = 3  #verbosity: errors(0), warnings(1), info(2), hints(3)
sc.logging.print_versions()



######### load andata #########
#path = '/data/zhanglab/kdunton/6samples_cluster/deepseq_3_clustering/DESC/cluster_by_genes/ourL23train/classifier-train_vs_control_3.8.22/cluster_TrCtrl_together_noDelta3.24.22/test_DESC_rerun/TrCtrl_noDelta_classifierGenesOnly.csv'
#path = '/data/zhanglab/kdunton/6samples_cluster/deepseq_3_clustering/DESC/cluster_by_genes/ourL23train/classifier-train_vs_control_3.8.22/cluster_TrCtrl_together_noDelta3.24.22/test_DESC_decimals/TrCtrl_noDelta_classifierGenesOnly_7dec.csv'
#path = '/data/zhanglab/kdunton/6samples_cluster/deepseq_3_clustering/DESC/cluster_by_genes/ourL23train/classifier-train_vs_control_3.8.22/cluster_TrCtrl_together_noDelta3.24.22/TrCtrl_noDelta_classifierGenesOnly.csv'

#print("Path:")
#print(path)
#adata = sc.read(path)
  #this reads a counts matrix into anndata object; has to have cells in rownames and genes in columns


######### preprocessing #########
#the counts in adata have already been normalized by highRankGenes.r

## view object
#print("adata initial")
#print(adata)

#adata.raw = adata


n_neighbors = [8,10,15,18,19,20,21,22,23,24,25]
#n_neighbors = [21,22,23,24,25]
n_neighbors = np.array(n_neighbors)
#louvain_resolution = [0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,
 #                     0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,
  #                    0.9]
louvain_resolution = [0.7,0.75,0.8,0.85,0.9]
louvain_resolution = np.array(louvain_resolution)

######### run DESC #########
for n in n_neighbors:
  print("n_neighbors=" + str(n))
  for res in louvain_resolution:
    print("louvain_resolution=" + str(res))

    #path = '/work/pi_yingzhang_uri_edu/kdunton/RNAseq/data/memonet_data/counts_matrices/L23_DGmtx_dec7.csv'
    path = '/work/pi_yingzhang_uri_edu/kdunton/RNAseq/cluster_by_genes/DESC_parameter_test/L23_0.3_DGmtx_dec7.csv'
    print("Path:"+ str(path))
    adata = sc.read(path)
    print("adata initial")
    print(adata)
    adata.raw = adata

    save_dir = ["result_DESC7.n"+str(n)+".L"+str(res)]
      #save_dir is a list
    #convert to string
    save_dir = ''.join(map(str,save_dir))

    adata = desc.train(adata, dims=[adata.shape[1], 64, 32], tol=0.001, n_neighbors=n,
                   batch_size=256, louvain_resolution=[res],
                   save_dir=save_dir, do_tsne=True, learning_rate=300,
                   do_umap=True, num_Cores_tsne=4,
                   save_encoder_weights=False)


    ######### Visualize #########
    #Umap: plot the clusterIDs
    basis = ["umap" + str(res)]
    basis = ''.join(map(str,basis))
    print(basis)
    color = ["desc_" + str(res)]
    color = ''.join(map(str,color))
    print(color)
    save = ["desc.dec7.n" + str(n) + ".L" + str(res) + ".png"]
    print(save)
    save = ''.join(map(str,save))

    sc.pl.scatter(adata, basis=basis, color=color, save = save)


######### save result #########
## save AnnData object into .h5ad file ##
#adata.write('desc_result.h5ad')

## save into csv ##
# tsne coordinates:
    #tsne = ["X_tsne" + str(res)]
    #tsne = ''.join(map(str,tsne))
    #obsm_data=pd.DataFrame(adata.obsm[tsne])
    #save_tsne = ["tsne_dec7.n" + str(n) + ".L" + str(res) + ".csv"]
    #save_tsne = ''.join(map(str,save_tsne))
    #obsm_data.to_csv(save_tsne,sep=",")

# umap coordinates:
    umap = ["X_umap" + str(res)]
    umap = ''.join(map(str,umap))
    obsm_data=pd.DataFrame(adata.obsm[umap])
    save_umap = ["umap_dec7.n" + str(n) + ".L" + str(res) + ".csv"]
    save_umap = ''.join(map(str,save_umap))
    obsm_data.to_csv(save_umap,sep=",")

# clustering result/labels:
    clu = ["desc_" + str(res)]
    clu = ''.join(map(str,clu))
    obs_data=pd.DataFrame(adata.obs[clu])
    save_clu = ["clusters_dec7.n" + str(n) + ".L" + str(res) + ".csv"]
    save_clu = ''.join(map(str,save_clu))
    obs_data.to_csv(save_clu,sep=",")


print("done!")
