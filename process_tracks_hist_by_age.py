"""
Consolidate Hood Canal dolphin DNA release experiments into a dictionary of snapshots from 
42 hc_dolph particle release experiments at a reference time, refT
"""
import netCDF4 as nc
import numpy as np
import pickle
from datetime import datetime, timedelta
import pandas as pd
import pytz

import sys
import os

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()

home = '/data2/pmr4/eab32/'
track_dir0 = home+'LO_output/tracks/'

D = dict()
D['k'] = 0.02/3600 #units: data 0.02 1/hr, multiply by hr/sec to get 1/sec

f_list = os.listdir(track_dir0)
f_list.sort()
f_list = [x for x in f_list if x[:11]=='hc_dolph_3d']


#save grid parameters
grid_fn = track_dir0+f_list[0]+'/grid.nc'
dsg = nc.Dataset(grid_fn)
D['lon']=dsg['lon_rho'][:]
D['lat']=dsg['lat_rho'][:]
D['mask']=dsg['mask_rho'][:]
D['h']=dsg['h'][:]


# establish histogram bins
# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

# MAP
# set domain limits
pad = .01
# plot region around delta pier, using release point as indicator
aa = [lon0-pad, lon0+pad,
   lat0-pad, lat0+pad]
nbins = 50
bin_lon_edges=np.linspace(aa[0], aa[1],nbins+1)
bin_lat_edges=np.linspace(aa[2], aa[3],nbins+1)
D['bin_lon_edges'] = bin_lon_edges
D['bin_lat_edges'] = bin_lat_edges
xx, yy = np.meshgrid(bin_lon_edges[:-1]+0.5*(bin_lon_edges[1]-bin_lon_edges[0]),bin_lat_edges[:-1]+0.5*(bin_lat_edges[1]-bin_lat_edges[0]))


# load eDNA concentrations
eDNA_fn = home+'LO_data/eDNA/ESP_Feb2023_hourly.csv'
df = pd.read_csv(eDNA_fn,sep=',',engine='python')
df.date = pd.to_datetime(df.date)

varlist = ['lon','lat','z'] # options: 'lon', 'lat', 'cs', 'ot', 'z', 'salt', 'temp', 'zeta', 'h', 'u', 'v', 'w'

rel_count=0
for f in f_list:
    
    track_dir = track_dir0+f
    
    # build a keyname from the release filename
    file_list = os.listdir(track_dir)
    file_list = [x for x in file_list if x[:3]=='rel']
    rel_fn = file_list[0]
    
    ds = nc.Dataset(track_dir+'/'+rel_fn)
    
    plon = ds['lon'][:].data
    plat = ds['lat'][:].data
    ot = ds['ot'][:].data
    nt = len(ot)
    
    if f==f_list[0]:
        nrel = len(f_list)
        
        nx,ny = xx.shape
        D['hist'] = np.zeros((nrel,nt,ny,nx))
        D['C0'] = np.zeros((nrel))
    
    
    T0 = datetime(1970,1,1,tzinfo=pytz.utc)+timedelta(seconds=ot[0]) #datetime of release
    D['C0'][rel_count] = df[df.date==T0].PB_quantity_mean.values[0] #look up the C0 for the release time
    
    
    for t in range(nt):
        decay = np.exp(-(ot[t]-ot[0])*D['k']) #scale by decay as a function of time in seconds since release
        hist = np.histogram2d(plon[t,:],plat[t,:],bins=[D['bin_lon_edges'],D['bin_lat_edges']])
        D['hist'][rel_count,t,:,:] = decay*D['C0'][rel_count]*hist[0].T
    
    ds.close()
    rel_count+=1

outfn = track_dir0 + 'all_3d_hc_dolph_hist_by_age.p'
pickle.dump(D,open(outfn,'wb'))