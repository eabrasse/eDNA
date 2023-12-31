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
import efun

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

plt.close('all')
home = '/data2/pmr4/eab32/'

atoz = string.ascii_lowercase

tab10 = plt.get_cmap('tab10',10)
Paired = plt.get_cmap('Paired',12)

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

# load background map
# figname = home+'LO_data/eDNA/bangor_gmaps.png'
# img_extent = (-122.754167, -122.695833, 47.718611, 47.788056)
axis = [-122.754167, -122.695833, 47.718611, 47.788056]

#find June sampling particle tracks
moor_fn = home+'LO_data/eDNA/Feb2023_moorings.p'

moor_dict = pickle.load(open(moor_fn,'rb'))
moor_list = moor_dict['moor_list']
dt_list = moor_dict['dt_list']
nmoor = len(moor_list)
rainbow = plt.get_cmap('rainbow',nmoor)
moor_lat_list = [moor['lat'] for moor in moor_list]
ns_inds = np.argsort(moor_lat_list)[::-1]

#get mean concentration 
data_fn = home+'LO_data/eDNA/ESP_Feb2023_hourly.csv'
df = pd.read_csv(data_fn,sep=',',engine='python')
DNA_mean = np.mean(df.PB_quantity_mean)

# with Image.open(figname) as img:

fig = plt.figure(figsize=(12, 8))
gs = GridSpec(nmoor*3,3)

ax_mean = plt.subplot(gs[:(nmoor),-1])
ax_var = plt.subplot(gs[(nmoor):(2*nmoor),-1])
ax_std2mean = plt.subplot(gs[(2*nmoor):,-1])

# ax = plt.subplot(gs[:,0],projection=ccrs.PlateCarree())
ax = plt.subplot(gs[:,0])

# # get grid from random file
grid_fn= home+ 'LO_output/tracks/all_3d_hc_dolph_releases.p'

D = pickle.load(open(grid_fn,'rb'))

# gather some fields, for convenience
lonr = D['metadata']['lon'][:]
latr = D['metadata']['lat'][:]
maskr = D['metadata']['mask'][:]


# MAP
xgrid,ygrid = efun.ll2xy(lonr,latr,lon0,lat0)

# ax.imshow(img, extent=img_extent,transform=ccrs.PlateCarree())
ax.pcolormesh(xgrid,ygrid,maskr,cmap=cmap_mask,shading='nearest',zorder=0)
ax.contour(xgrid,ygrid,maskr,levels=[0.5],colors=['k'],linewidths=[1.5])
ax.set_xlabel('Dist from pen (m)')
ax.set_ylabel('Dist from pen (m)')
# gl=ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)
# gl.top_labels = False
# gl.right_labels = False

moor_count=0
mooring_axes = []

for ind in ns_inds:

    moor = moor_list[ind]
    
    x,y = efun.ll2xy(moor['lon'],moor['lat'],lon0,lat0)
    ax.plot(x,y,marker='d',mec='k',mfc = rainbow(moor_count),markersize=10)

    axm = plt.subplot(gs[(3*moor_count):(3*(moor_count+1)),1])

    # moor['const_particle_bin'][moor['const_particle_bin']==0] = np.nan
    axm.plot(dt_list,DNA_mean*moor['const_particle_bin'],linestyle='dashed',color=rainbow(moor_count))
    # moor['TV_particle_bin'][moor['TV_particle_bin']==0] = np.nan
    axm.plot(dt_list,moor['TV_particle_bin'],linestyle='solid',color=rainbow(moor_count))

    axm.set_yscale('log')
    axm.text(0.1,0.9,'{}) {}'.format(atoz[moor_count],moor['label']),color='k',transform=axm.transAxes,ha='left',va='top')
    axm.grid()

    # get universal ylims
    if moor_count==0:
        ymax = max([np.max(DNA_mean*moor['const_particle_bin']),np.max(moor['TV_particle_bin'])])
        ymin = min([np.min(DNA_mean*moor['const_particle_bin']),np.min(moor['TV_particle_bin'])])
    ymax = max([np.max(DNA_mean*moor['const_particle_bin']),np.max(moor['TV_particle_bin']),ymax])
    ymin = min([np.min(DNA_mean*moor['const_particle_bin']),np.min(moor['TV_particle_bin']),ymin])

    # format date/time axis
    if moor_count<(nmoor-1):
        axm.set_xlabel('')
        axm.set_xticklabels(['' for xtl in axm.get_xticklabels()])
    else:
        axm.set_xlabel('Date & time')
        axm.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d %H:%m"))
        plt.setp( axm.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')

    mooring_axes.append(axm)

    # calculate stats
    
    moor['dist_from_pen'] = np.sqrt(x**2+y**2)
    const_std = np.std(DNA_mean*moor['const_particle_bin'][moor['const_particle_bin']>0])
    const_mean =np.mean(DNA_mean*moor['const_particle_bin'][moor['const_particle_bin']>0])
    moor['const_mean'] = const_mean
    moor['const_var'] = const_std**2
    moor['const_std2mean'] = const_std/const_mean
    TV_mean = np.mean(moor['TV_particle_bin'][moor['TV_particle_bin']>0])
    TV_std = np.std(moor['TV_particle_bin'][moor['TV_particle_bin']>0])
    moor['TV_mean'] = TV_mean
    moor['TV_var'] = TV_std**2
    moor['TV_std2mean'] = TV_std/TV_mean

    ax_mean.plot(moor['dist_from_pen'],moor['const_mean'],marker='o',linestyle='None',mfc='none',mec=rainbow(moor_count))
    ax_mean.plot(moor['dist_from_pen'],moor['TV_mean'],marker='o',linestyle='None',mfc=rainbow(moor_count),mec='None')
    ax_var.plot(moor['dist_from_pen'],moor['const_var'],marker='o',linestyle='None',mfc='none',mec=rainbow(moor_count))
    ax_var.plot(moor['dist_from_pen'],moor['TV_var'],marker='o',linestyle='None',mfc=rainbow(moor_count),mec='None')
    ax_std2mean.plot(moor['dist_from_pen'],moor['const_std2mean'],marker='o',linestyle='None',mfc='none',mec=rainbow(moor_count))
    ax_std2mean.plot(moor['dist_from_pen'],moor['TV_std2mean'],marker='o',linestyle='None',mfc=rainbow(moor_count),mec='None')

    moor_count+=1


# ax.axis([-122.735,-122.725,47.74,47.75])
lonaxes = np.array([-122.735,-122.725]);lataxes=np.array([47.74,47.75])
xaxes,yaxes = efun.ll2xy(lonaxes,lataxes,lon0,lat0)
ax.axis([xaxes[0],xaxes[1],yaxes[0],yaxes[1]])
ax.set_aspect(1)

ax_mean.text(0.1,0.9,'{}) Mean of nonzeros'.format(atoz[moor_count]),color='k',transform=ax_mean.transAxes,ha='left',va='top')
ax_var.text(0.1,0.9,'{}) Variance of nonzeros'.format(atoz[moor_count+1]),color='k',transform=ax_var.transAxes,ha='left',va='top')
ax_std2mean.text(0.1,0.9,'{}) Ratio of standard dev to mean'.format(atoz[moor_count+2]),color='k',transform=ax_std2mean.transAxes,ha='left',va='top')

for ax_stat in ax_mean,ax_var:
    ax_stat.set_xlabel('')
    ax_stat.set_xticklabels(['' for xtl in ax_stat.get_xticklabels()])
    ax_stat.set_yscale('log')
    ax_stat.set_ylabel('Part wt')

ax_std2mean.set_xlabel('Dist from pen (m)')
for ax_stat in ax_mean,ax_var,ax_std2mean:
    ax_stat.grid()
    

for axm in mooring_axes:
    axm.set_ylim([ymin,ymax])
    axm.set_ylabel('Part wt')

mooring_axes[0].text(0.9,0.9,'solid=time varying\ndashed=constant',transform=mooring_axes[0].transAxes,ha='right',va='top',fontsize=8)
ax_mean.text(0.95,0.8,'filled=time varying\nempty=constant',transform=ax_mean.transAxes,ha='right',va='top',fontsize=8)

fig.subplots_adjust(right=0.98,left=0.1,bottom=0.1,top = 0.98)
outfn = home+'etools/plots/Feb2023_mooringsF.png'
plt.savefig(outfn)
print(f'Saved to {outfn}')
plt.close()


