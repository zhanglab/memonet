# this script was written by Dr. EunJung Hwang
# it runs a smoothing function on the silhouette scores of each parameter combo to determine the peak silhouette score, and thus, the settings to use for DESC clustering

# to use with Jupyter notebook, run 'pip install notebook' in terminal, then 'jupyter notebook' to run line-by-line

import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
from pathlib import Path


SS=pd.read_csv('~/Downloads/RNAseq/cluster_by_genes/DESC_parameter_test/parameter_silhouette_score_avgs.csv')
param1=[]
param2=[]
sp1=[]
sp2=[]
for i in range(len(SS)):
    params = SS.parameters.values[i]
    ind = params.find('.')
    param1.append(float(params[1:ind]))
    param2.append(float(params[ind+2:]))
    sp1.append(params[1:ind])
    sp2.append(params[ind+2:])
uP1 = np.unique(param1)
uP2 = np.unique(param2)
ss_list =[]
ss_mat = np.ones((len(uP1),len(uP2)))
for i,u1 in enumerate(uP1):
    for j, u2 in enumerate(uP2):
        ind1 = np.flatnonzero(param1==u1)
        ind2 = np.flatnonzero(param2==u2)
        ind = np.intersect1d(ind1,ind2)
        ss_mat[i,j] = np.mean(SS.overall_silhouette_score.values[ind])


img = ndimage.gaussian_filter(ss_mat, sigma=(2, 1), order=0)
inds=np.unravel_index(np.argmax(img, axis=None), img.shape)

# left plot (regular heatmap):
fig, axes = plt.subplots(1,2,figsize=[10,5])
plt.subplot(1,2,1)
plt.matshow(ss_mat, fignum=0)
plt.plot(inds[1],inds[0],color='red',marker='x')
plt.xticks(np.arange(len(uP2)),labels=np.unique(sp2))
plt.xlabel('Louvain resolution')
plt.ylabel('N_neighbor')
plt.yticks(np.arange(len(uP1)),labels=np.unique(sp1))
plt.gca().invert_yaxis()
plt.subplot(1,2,2)

# right plot (smooth):
print(uP1[inds[0]], uP2[inds[1]])
plt.matshow(img, fignum=0)
plt.plot(inds[1],inds[0],color='red',marker='x')
plt.xticks(np.arange(len(uP2)),labels=np.unique(sp2))
plt.xlabel('Louvain resolution')
plt.ylabel('N_neighbor')
plt.yticks(np.arange(len(uP1)),labels=np.unique(sp1))
plt.gca().invert_yaxis()

# save plots:
plt.savefig('Smooth_Silhouette.png')
plt.savefig('Smooth_Silhouette.svg')
