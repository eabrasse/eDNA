import os
import sys
import numpy as np
from datetime import datetime, timedelta
import pytz
import netCDF4 as nc
import pandas as pd
import efun
import pickle
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmocean as cmo
from matplotlib.gridspec import GridSpec

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)

tab10 = plt.get_cmap('tab10',10)

home = '/data2/pmr4/eab32/'

#find June sampling particle tracks
track_dir0 = home+'LO_output/tracks/'

f_list = os.listdir(track_dir0)
f_list.sort()
f_list = [x for x in f_list if (x[:11]=='hc_dolph_3d')&(x[-4]=='6')]

#arbirtrary extent of interest, based on June 5 sampling
aa = [-122.754167, -122.695833, 47.718611, 47.788056]
gridfn = track_dir0+f_list[0]+'/grid.nc'
dsg = nc.Dataset(gridfn)
lonr = dsg['lon_rho'][:]
latr = dsg['lat_rho'][:]
lonvec = lonr[0,:]
latvec = latr[:,0]
nlon = np.shape(lonvec[(lonvec>aa[0])&(lonvec<aa[1])])[0]
nlat = np.shape(latvec[(latvec>aa[2])&(latvec<aa[3])])[0]
# nbins = 100
bin_lon_edges=np.linspace(aa[0], aa[1], nlon+1)
bin_lat_edges=np.linspace(aa[2], aa[3], nlat+1)
xx, yy = np.meshgrid(bin_lon_edges[:-1]+0.5*(bin_lon_edges[1]-bin_lon_edges[0]),bin_lat_edges[:-1]+0.5*(bin_lat_edges[1]-bin_lat_edges[0]))
x,y = efun.ll2xy(xx,yy,lonr[0,0],latr[0,0])
dx = np.mean(np.diff(x,axis=1))
dy = np.mean(np.diff(y,axis=0))

#set reference time
Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

first_sample_utc = Pacific.localize(datetime(2023,6,5,10,30)).astimezone(utc)
first_release = first_sample_utc-timedelta(days=2)
last_sample_utc = Pacific.localize(datetime(2023,6,5,15,30)).astimezone(utc)
dt_list0 = pd.date_range(start=first_sample_utc,end = last_sample_utc, freq="15min").to_pydatetime().tolist()

# radius = 100
depth = 5

dt = Pacific.localize(datetime(2023,6,5,11,0)).astimezone(utc)
dtts = pd.Timestamp(dt) #useful for comparisons with pandas dataframes

concentration = np.zeros(xx.shape)

nfile = len(f_list)
filecount=1
# unlike plotting code, open netCDF files only once and loop through timesteps within files
for f in f_list:

    track_dir = track_dir0+f
    
    print(f'Working on file {filecount} of {nfile}...')

    # build a keyname from the release filename
    file_list = os.listdir(track_dir)
    file_list = [x for x in file_list if x[:3]=='rel']
    rel_fn = file_list[0]

    ds = nc.Dataset(track_dir+'/'+rel_fn)

    ot = ds['ot'][:].data
    
    dt_list = [pd.Timestamp(utc.localize(datetime(1970,1,1)+timedelta(seconds=ott))) for ott in ot]

    ot0 = dt_list[0]
    ot1 = dt_list[-1]
    
    if (dt>ot0)&(dt<ot1):
        t = argnearest(dt_list,dtts)
        zmask = ds['z'][t,:]>(ds['zeta'][t,:]-depth)
        if np.sum(zmask)>0:
            # assume there are fewer particles > 5m deep than grid cells
            # 
            # to find the particle concentration within 100m and above 5m depth centered at each grid cell,
            # add each particle to the grid cell centers within 100m 
            lonp = ds['lon'][t,zmask]
            latp = ds['lat'][t,zmask]
            
            hist = np.histogram2d(lonp,latp,bins=[bin_lon_edges,bin_lat_edges])
            concentration+=hist[0].T
           
            # for particle in range(len(lonp)):
            #     x,y = efun.ll2xy(lonr,latr,lonp[particle],latp[particle])
            #     lldist = np.sqrt(x**2+y**2)
            #     lldistma = lldist<radius
            #     concentration+=lldistma
            
    ds.close()
    filecount+=1
volume = depth*dx*dy
concentration = concentration/volume

fig = plt.figure(figsize=(12, 8))
gs = GridSpec(2,3)
ax = fig.add_subplot(gs[:,:2])

# colorscale plotting parameters
# dyemax = 5e-2
# dyemin = 1e-4
# vmax = dyemax
# vmin = dyemin
colmap = cmo.cm.matter
# normal = matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax)
normal = matplotlib.colors.LogNorm()

p=ax.pcolormesh(xx,yy,concentration,cmap = colmap,norm=normal) # plot concentration
maskr = dsg['mask_rho'][:]
ax.contour(lonr,latr,maskr,levels=[0.5],colors='k',linewidths=1,linestyles='solid') # plot shoreline
cbaxes = inset_axes(ax, width="6%", height="40%", loc='lower left',bbox_transform=ax.transAxes,bbox_to_anchor=(0.,0,1,1))
ax.axis(aa)
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical') # add concentration colorbar

axlon = fig.add_subplot(gs[0,2])
lon_conc = np.sum(concentration,axis=0)
axlon.bar(xx[0,:],lon_conc,color=tab10(0),width=0.001)
axlon.set_xlabel('Longitude')
axlon.set_ylabel('Sum particle concentration')

axlat = fig.add_subplot(gs[1,2])
lat_conc = np.sum(concentration,axis=1)
axlat.bar(yy[:,0],lat_conc,color=tab10(1),width=0.001)
axlat.set_xlabel('Latitude')
axlat.set_ylabel('Sum particle concentration')

outfn = home+'etools/plots/HCJune5_concentration_binned.png'
plt.savefig(outfn)
print(f'Saved to {outfn}')
plt.close()