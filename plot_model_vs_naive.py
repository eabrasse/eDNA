#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot dolphin pen particles nearby June 5 sampling stations
"""

import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
import efun
import cmocean
import pickle

# some helpful plotting commands
plt.close('all')


home = '/Users/elizabethbrasseale/Projects/eDNA/'
data_fn = home+'data/HCJune5_sampling_station_neighbors.p'
D = pickle.load(open(data_fn,'rb'))
rainbow = plt.get_cmap('rainbow',len(D.keys()))

sampling_data_fn = home+'data/June5_sampling_locations.csv'
dfs = pd.read_csv(sampling_data_fn,sep=',',engine='python')

radius =100
depth = 5
overvol = 1/(depth*np.pi*radius**2)
# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773
model_conc = np.zeros(len(D.keys()))
eDNA_present = np.zeros(len(D.keys()))
dist = np.zeros(len(D.keys()))
figure,axs = plt.subplots(nrows=1,ncols=2,figsize=(9,4))
count=0
for fieldID in D.keys():
    if count==0:
        radius_list = D[fieldID]['radius_list']
        rind = radius_list.index(radius)
        depth_list = D[fieldID]['depth_list']
        dind = depth_list.index(depth)
    g = dfs[dfs.sampling_location==fieldID]
    
    model_conc[count] = D[fieldID]['station_neighbors'][dind,rind]*overvol
    eDNA_present[count] = g['eDNA_present']
    # axs[0].plot(model_conc,g['eDNA_present'],linestyle='none',marker='o',mfc=rainbow(count),mec=rainbow(count))
    
    x,y = efun.ll2xy(g['long'].values[0],g['lat'].values[0],lon0,lat0)
    dist[count]=np.sqrt(x**2+y**2)
    # axs[1].plot(dist,g['eDNA_present'],linestyle='none',marker='o',mfc=rainbow(count),mec=rainbow(count))
    
    count+=1
nbins = 5
for ax,var in [axs[0],model_conc],[axs[1],dist]:
    bin_edges = np.linspace(var.min(),var.max(),nbins+1)
    dbin = 0.5*(bin_edges[1]-bin_edges[0])
    bins = dbin+bin_edges[:-1]
    fract_in_bin = np.zeros(nbins)
    for i in range(nbins):
        var_inds = (var>bin_edges[i])&(var<bin_edges[i+1])
        fract_in_bin[i] = np.mean(eDNA_present[var_inds])
    ax.bar(bins,fract_in_bin,width=dbin)
    ax.grid()
    ax.set_xticks(bins)

axs[0].set_xlabel('Model concentration')
axs[0].set_ylabel('eDNA fraction')
axs[0].set_yticks([0,0.25,0.5,0.75,1])
# axs[0].set_yticklabels(['Absent','Present'])
axs[0].text(0.1,0.7,'a) Model particle concentration\nvs sampled eDNA detections',transform=axs[0].transAxes)

axs[1].set_xlabel('Dist (m)')
axs[1].set_ylabel('')
axs[1].set_yticks([0,0.25,0.5,0.75,1])
axs
axs[1].set_yticklabels(['','','','',''])
axs[1].text(0.1,0.7,'b) Naïve distance\nvs sampled eDNA detections',transform=axs[1].transAxes)

plt.show(block=False)
plt.pause(0.1)
# outfn = home+'/plots/model_vs_naive.png'
# plt.savefig(outfn)
# print(f'Saved to {outfn}')