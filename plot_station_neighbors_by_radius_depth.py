#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot dolphin pen particles nearby June 5 sampling stations
"""

import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from PIL import Image
import numpy as np
from datetime import datetime, timedelta
import pytz
import netCDF4 as nc
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.dates as mdates
import efun
import math
import cmocean
import pickle

# some helpful plotting commands
plt.close('all')
home = '/data2/pmr4/eab32/'
data_fn = home+'LO_data/eDNA/HCJune5_sampling_station_neighbors.p'
D = pickle.load(open(data_fn,'rb'))

nstation = len(D.keys())
ncols = 5
nrows = math.floor(nstation/ncols)

fig,axs = plt.subplots(nrows,ncols,figsize=(12, 8))

station_count = 0
ylim_list = []
noyes_colorlist = ['red','green']
for fieldID in D.keys():
    
    ax = axs.flatten()[station_count]
    
    station = D[fieldID]
    station['station_neighbors'][station['station_neighbors']==0]+=1
    
    if station_count==0:
        nradius = len(station['radius_list'])
        ndepth = len(station['depth_list'])
        cmo_deep = plt.get_cmap('inferno', ndepth*2)
        ax.text(0.1,0.8,'Depth threshold',transform=ax.transAxes,color='k',ha='left',fontweight='normal')
        for d in range(ndepth):
            yloc = 0.7-0.1*d
            ax.text(0.1,yloc,'{} m'.format(station['depth_list'][d]),color=cmo_deep(d),transform=ax.transAxes,ha='left',fontweight='bold')
    
    for d in range(ndepth):
        ax.loglog(station['radius_list'],station['station_neighbors'][d,:],mfc=cmo_deep(d),linestyle='None',marker='o',mec='none')
    
    ax.text(0.1,0.9,f'Station {fieldID}',transform=ax.transAxes,color=noyes_colorlist[station['eDNA_present']])
    ylim_list.append(ax.get_ylim()[1])
    
    station_count+=1

for ax in axs[:,0]:
    ax.set_ylabel('Particle count')
for ax in axs[-1,:]:
    ax.set_xlabel('Dist threshold (m)')

yl0 = 0.5
yl1 = max(ylim_list)
for ax in axs.flatten():
    
    ax.set_ylim([yl0,yl1])
    # ax.set_yscale('log')
    
outfn = home+'etools/plots/June5_sampling_station_neighbors.png'
plt.savefig(outfn)
print(f'Saved to {outfn}')