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
import pickle
import string

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)

plt.close('all')
home = '/data2/pmr4/eab32/'

atoz = string.lower_ascii

tab10 = plt.get_cmap('tab10',10)
Paired = plt.get_cmap('Paired',12)

# load background map
figname = home+'LO_data/eDNA/bangor_gmaps.png'
img_extent = (-122.754167, -122.695833, 47.718611, 47.788056)

#find June sampling particle tracks
moor_fn = home+'LO_data/eDNA/Feb2023_moorings.p'

moor = pickle.load(open(moor_fn,'rb'))
moor_list = moor['moor_list']
nmoor = len(moor_list)
rainbow = plt.get_cmap('rainbow',nmoor)
moor_lat_list = [moor['lat'] for moor in moor_list]
ns_inds = moor_lat_list.argsort()[::-1]

with Image.open(figname) as img:

    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(nmoor,3)

    ax = plt.subplot(gs[:,0],projection=ccrs.PlateCarree())

    ax.imshow(img, extent=img_extent,transform=ccrs.PlateCarree())
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    gl=ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    
    moor_count=0
    mooring_axes = []
    
    for ind in ns_inds:

        moor = moor_list[ind]
        
        ax.plot(moor['lon'],moor['lat'],marker='d',mec='k',mfc = rainbow(moor_count),markersize=12)
        
        axm = plt.subplot(gs[moor_count,1:])
        axm.plot(moor['dt_list'],moor['const_particle_bin'],linestyle='dashed',color=rainbow(moor_count))
        axm.plot(moor['dt_list'],moor['TV_particle_bin'],linestyle='solid',color=rainbow(moor_count))
        
        axm.set_yscale('log')
        axm.text(0.1,0.9,'{}) {}'.format(atoz[moor_count],moor['label']),color=rainbow(moor_count),transform=axm.transAxes,ha='left',va='top')
        axm.grid()
        
        # get universal ylims
        if moor_count==0:
            ymax = max([np.max(moor['const_particle_bin']),np.max(moor['TV_particle_bin'])])
            ymin = min([np.min(moor['const_particle_bin']),np.min(moor['TV_particle_bin'])])
        ymax = max([np.max(moor['const_particle_bin']),np.max(moor['TV_particle_bin']),ymax])
        ymin = min([np.min(moor['const_particle_bin']),np.min(moor['TV_particle_bin']),ymin])
        
        # format date/time axis
        if moor_count==(nmoor-1):
            axm.set_xlabel('')
            axm.set_xticklabels(['' for xtl in axm.get_xticklabels()])
        else:
            axm.set_xlabel('Date & time')
            axm.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d %H:%m"))
            plt.setp( axm.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
        
        mooring_axes.append(axm)
        moor_count+=1
        

    ax.axis(img_extent)

    for axm in mooring_axes:
        axm.set_ylim([ymin,ymax])

    outfn = home+'etools/plots/Feb2023_moorings.png'
    plt.savefig(outfn)
    print(f'Saved to {outfn}')
    plt.close()


