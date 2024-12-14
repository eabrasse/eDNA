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
from matplotlib import ticker
import matplotlib.dates as mdates
import pickle
import string
import efun
import pytz
from scipy.stats import linregress
import cmocean as cmo
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

tz = pytz.utc
pdt = pytz.timezone('America/Vancouver')

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('#526c30ff')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

blue_circle_color = '#37abc8ff'

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'

tab10 = plt.get_cmap('tab10',10)
Paired = plt.get_cmap('Paired',12)


#find June sampling particle tracks
data_fn = home+'data/Feb2023_moor_radius_counts.p'

D = pickle.load(open(data_fn,'rb'))
#D.keys() = ['moor_list','particle_age_bins','radius_list','sample_dt','zref']
for key in D.keys():
    locals()[key] = D[key]

volume_list = [(np.pi*rad**2)*np.abs(zref) for rad in radius_list]
over_volume_list = [1/vol for vol in volume_list]
nrad = len(radius_list)
nmoor = len(moor_list)

agescale = 1/3600
for m in range(nmoor):
    
    moor = moor_list[m]
    
    fw,fh = efun.gen_plot_props()
    fig = plt.figure(figsize=(fw*2,fh))
    gs = GridSpec(1,2)
    ax_conc = plt.subplot(gs[0])
    ax_age = plt.subplot(gs[1])
    
    conc_list = [moor['count'][r]*over_volume_list[r] for r in range(nrad)]
    
    ax_conc.plot(radius_list,conc_list,marker='o',mfc=blue_circle_color,alpha=1,mec='None')
    ax_age.plot(radius_list,particle_age_bins[m,:],marker='o',mfc=tab10(1),alpha=1,mec='None')
    
    ax_conc.text(0.1,0.9,moor['label'],transform=ax_conc.transAxes)
    ax_conc.set_xlabel('radius (m)')
    ax_conc.set_ylabel('concentration')
    
    ax_age.set_xlabel('radius (m)')
    ax_age.set_ylabel('mean age')
    

    fig.subplots_adjust(right=0.98,left=0.15,bottom=0.15,top = 0.98,wspace=0.3)

    plt.show(block=False)
plt.pause(0.1)
# plt.close()


