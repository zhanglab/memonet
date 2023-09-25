# this script overlays individual mouse proportions as dots on top of the bar graph of train/control proportions across clusters

import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
from pathlib import Path

SS = pd.read_csv(‘~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/figures/ClusterMouseProp_mouseDots.csv’)
plt.rcParams[‘svg.fonttype’] = ‘none’ # do not convert text to path
plt.rcParams[‘font.sans-serif’] = ‘Arial’
plt.rcParams[‘font.size’] = 6
plt.rcParams[‘font.stretch’] = ‘normal’
figure = plt.figure(figsize=np.array([2.2,2.5]))
axe = plt.axes()
axe.spines[‘top’].set_visible(False)
axe.spines[‘right’].set_visible(False)
colors = [‘Black’, ‘DarkGreen’,‘darkGray’,‘Green’,‘LightGray’,‘LimeGreen’]
colors2 = [‘dimgray’,‘lime’]
for stim in np.array([0,1]):
    stack_val = []
    for mouse in np.arange(stim+1,stim+6,2):
        prop = SS.mouseProp[np.flatnonzero(SS.mouse == mouse)].values
        stack_val.append(prop)
        plt.plot(np.arange(6)*2 + 0.5*stim-0.25, prop, ‘o’,color=colors[mouse-1], markersize=3.2)
    plt.bar(np.arange(6)*2 + 0.5*stim-0.25, np.mean(stack_val,axis=0), color=colors2[stim], width=0.5, edgecolor=‘black’, linewidth=0.5)
plt.xticks(2*np.arange(6),labels=[‘C0’,‘C1’,‘C2’,‘C3’,‘C4’,‘C5’])
plt.xlabel(‘Cluster’, fontsize=8)
plt.ylabel(‘Proportion’, fontsize=8)
plt.savefig(‘~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/figures/ClusterMouseProp_mouseDots.svg’)


