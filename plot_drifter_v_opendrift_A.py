import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
# from datetime import datetime, timedelta
# import pytz
from matplotlib.gridspec import GridSpec
# import matplotlib.dates as mdates
import pickle
from matplotlib import ticker
import string
import efun
import cmocean as cmo
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
from datetime import datetime,timedelta
import netCDF4 as nc
import pytz
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
from PIL import Image
from matplotlib.collections import LineCollection

pdt = pytz.timezone('America/Vancouver')
utc = pytz.utc


# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])
tab10 = plt.get_cmap('tab10',10)
Set2 = plt.get_cmap('Set2',8)
atoz = string.ascii_lowercase
rainbow = plt.get_cmap('rainbow')

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'

# load background map
figname = home+'data/bangor_satellite_map.png'
img_extent = (-122.743707, -122.726679, 47.736877, 47.744494)
a_ll = (-122.740707, -122.726679, 47.736877, 47.744494)

#dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

sample_day = datetime(2024,10,16)
# load in sampling locations
drifter_fn = home+'data/MSdrifter_data_10.16.2024.proc.txt'
drifter_df = pd.read_csv(drifter_fn,sep=',',engine='python')
inds = range(len(drifter_df.index))
drifter_df['datetime'] = [datetime(drifter_df['Year'].values[i],drifter_df['    Month'].values[i],drifter_df['    Day'].values[i])+timedelta(seconds=int(drifter_df['    Seconds past Day'].values[i])) for i in inds]
drifter_df['datetime'] = drifter_df['datetime'].dt.tz_localize('utc').dt.tz_convert('America/Vancouver')


meta_fn = home+'data/Drifter metadata.csv'
meta_df = pd.read_csv(meta_fn,sep=',',engine='python')
meta_df['dt0'] = [datetime(2024,10,16,int(meta_df['Release time (PDT)'].values[i][:2]),int(meta_df['Release time (PDT)'].values[i][3:])) for i in range(len(meta_df.index))]
meta_df['dt0'] = meta_df['dt0'].dt.tz_localize('America/Vancouver')
meta_df['dt1'] = [datetime(2024,10,16,int(meta_df['Recovery time (PDT)'].values[i][:2]),int(meta_df['Recovery time (PDT)'].values[i][3:])) for i in range(len(meta_df.index))]
meta_df['dt1'] = meta_df['dt1'].dt.tz_localize('America/Vancouver')
# meta_df['Release x'],meta_df['Release y'] = efun.ll2xy(meta_df['Release Lon'],meta_df['Release Lat'],lon0,lat0)
meta_df['Recovery x'],meta_df['Recovery y'] = efun.ll2xy(meta_df['Recovery Lon'],meta_df['Recovery Lat'],lon0,lat0)

relnum = 5
lag = 0
opendrift_fn = home+'data/opendrift/HCdrifter_DPmask_rel{}_{}hrlag_verticalmotion_NOAAwind.nc'.format(relnum,lag)
ds = nc.Dataset(opendrift_fn)
nseed,ntime = ds['lon'].shape
pt = ds['time'][:]
pt0 = utc.localize(datetime(1970,1,1)+timedelta(seconds=pt[0])).astimezone(pdt)
pt1 = utc.localize(datetime(1970,1,1)+timedelta(seconds=pt[-1])).astimezone(pdt)
px1, py1 = efun.ll2xy(ds['lon'][:,-1],ds['lat'][:,-1],lon0,lat0)


# for rel in meta_df['Release #']:
rel = meta_df['Release #'][relnum]
dt0 = meta_df['dt0'].loc[rel]
dt1 = meta_df['dt1'].loc[rel]
driftx1 = meta_df['Recovery x'].loc[rel]
drifty1 = meta_df['Recovery y'].loc[rel]
pr1 = np.sqrt((px1-driftx1)**2+(py1-drifty1)**2) # distance of particle end points from drifter end points

with Image.open(figname) as img:

    fw,fh = efun.gen_plot_props()
    fig = plt.figure(figsize=(fw*2.5,fh*1.))
    gs = GridSpec(1,1)

    ax = plt.subplot(gs[0,0],projection=ccrs.PlateCarree())
    ax.imshow(img, extent=img_extent,transform=ccrs.PlateCarree())
    # gl=ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)
    # gl.top_labels = False
    # gl.right_labels = False
    
    # #open drift prediction
    for i in range(nseed):
        points = np.array([ds['lon'][i,:], ds['lat'][i,:]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments,linestyle='solid',color='w',alpha=1,linewidth=0.1,zorder=100)
        line = ax.add_collection(lc)
    # ax.plot(ds['lon'][0,0],ds['lat'][0,0],marker='x',lw=3,markersize=8,mec='w',mfc='w',zorder=5000)
    # ax.scatter(ds['lon'][:,0],ds['lat'][:,0],zorder=50,s=10,alpha=0.3,edgecolors='w',c='none',lw=0.5,)
    # latref= 47.741
    rmask = pr1<100
    rinds = np.argwhere(rmask)
    for ind in rinds:
        points = np.array([ds['lon'][ind[0],:], ds['lat'][ind[0],:]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments,linestyle='dashed',color=Set2(int(rel)),alpha=1,linewidth=1,zorder=100)
        line = ax.add_collection(lc)
        ax.plot(ds['lon'][ind,0],ds['lat'][ind,0],marker='o',lw=3,markersize=5,mec=Set2(int(rel)),mfc='w',zorder=5000)
    
    #drifter tracks
    ddf = drifter_df[(drifter_df['datetime']>dt0)&(drifter_df['datetime']<dt1)]
    
    lats = np.array(ddf['   Latitude'].values[:])
    if rel>0:
        lats = np.insert(lats,0,meta_df['Release Lat'].loc[rel])
    lats = np.append(lats,meta_df['Recovery Lat'].loc[rel])
    
    lons = np.array(ddf['  Longitude'].values[:])
    if rel>0:
        lons = np.insert(lons,0,meta_df['Release Lon'].loc[rel])
    lons = np.append(lons,meta_df['Recovery Lon'].loc[rel])
    
    if rel>0:
        ax.plot(meta_df['Release Lon'].loc[rel],meta_df['Release Lat'].loc[rel],linestyle='none',marker='v',mec='k',mfc=Set2(int(rel)),markersize=10,zorder=350)
    ax.plot(lons,lats,linestyle='solid',marker='none',color=Set2(int(rel)),lw=1.5,zorder=350)
    ax.plot(ddf['  Longitude'],ddf['   Latitude'],linestyle='none',marker='o',mec=Set2(int(rel)),color=Set2(int(rel)),mfc=Set2(int(rel)),markersize=5,zorder=330)
    ax.plot(meta_df['Recovery Lon'].loc[rel],meta_df['Recovery Lat'].loc[rel],linestyle='none',marker='^',mec='k',mfc=Set2(int(rel)),markersize=10,zorder=360)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    ax.axis(a_ll)
    
    textstr0 = 'Drifter Release #{}\nStation: {}\nTime in: {}\nTime out: {}'.format(rel,meta_df['Station'].loc[rel],dt0.strftime("%H:%M PDT"),dt1.strftime("%H:%M PDT"))
    ax.text(0.05,0.95,textstr0,transform=ax.transAxes,ha='left',va='top',color=Set2(int(rel)))
    textstr1 = 'Particles\n\nTime in: {}\nTime out: {}'.format(pt0.strftime("%H:%M PDT"),pt1.strftime("%H:%M PDT"))
    ax.text(0.3,0.95,textstr1,transform=ax.transAxes,ha='left',va='top',color='w')

    # fig.subplots_adjust(right=0.98,left=0.15,bottom=0.11,top = 0.98)
    fig.tight_layout()
    # outfn = home+f'plots/drifter2024.10.16/drifter_GPS_R{rel:}_vs_opendrift.png'
    # plt.savefig(outfn)
    plt.show(block=False)
    plt.pause(0.1)
    # plt.close()


