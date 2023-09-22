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
D['metadata'] = dict()
D['metadata']['refT'] = datetime(2023,2,2,0,0,0,tzinfo=pytz.utc) #time stamp in UTC for experiment consolidation
D['metadata']['k'] = 0.02/3600 #units: data 0.02 1/hr, multiply by hr/sec to get 1/sec

f_list = os.listdir(track_dir0)
f_list.sort()
f_list = [x for x in f_list if x[:11]=='hc_dolph_3d']

#save grid parameters
grid_fn = track_dir0+f_list[0]+'/grid.nc'
dsg = nc.Dataset(grid_fn)
D['metadata']['lon']=dsg['lon_rho'][:]
D['metadata']['lat']=dsg['lat_rho'][:]
D['metadata']['mask']=dsg['mask_rho'][:]
D['metadata']['h']=dsg['h'][:]

# load eDNA concentrations
eDNA_fn = home+'LO_data/eDNA/ESP_Feb2023_hourly.csv'
df = pd.read_csv(eDNA_fn,sep=',',engine='python')
df.date = pd.to_datetime(df.date)

varlist = ['lon','lat','z'] # options: 'lon', 'lat', 'cs', 'ot', 'z', 'salt', 'temp', 'zeta', 'h', 'u', 'v', 'w'

for f in f_list:
    track_dir = track_dir0+f
    file_list = os.listdir(track_dir)
    file_list = [x for x in file_list if x[:3]=='rel']
    rel_fn = file_list[0]
    rel_date_str = f[f.index('2023'):]
    rel_hour_str = f[f.index('sh'):f.index('_2023')]
    keyname = 'rel_'+rel_date_str+'_'+rel_hour_str
    D[keyname] = dict()
    
    ds = nc.Dataset(track_dir+'/'+rel_fn)
    
    ot = ds['ot'][:].data
    D[keyname]['T0'] = datetime(1970,1,1,tzinfo=pytz.utc)+timedelta(seconds=ot[0]) #datetime of release
    td_list = [np.abs(datetime(1970,1,1,tzinfo=pytz.utc)+timedelta(seconds=t)-D['metadata']['refT']) for t in ot]
    refT_ind = td_list.index(min(td_list))
    
    # look at values at the index of the reference time
    D[keyname]['deltaT'] = ot[refT_ind]-ot[0] # number of seconds between release and refT
    for var in varlist:
        D[keyname][var] = ds[var][refT_ind,:]
    D[keyname]['C0'] = df[df.date==D[keyname]['T0']].PB_quantity_mean.values[0] #look up initial concentration near dolphin pen at release time
    D[keyname]['decay'] = np.exp(-D['metadata']['k']*D[keyname]['deltaT'])
    
    ds.close()

outfn = track_dir0 + 'all_3d_hc_dolph_releases.p'
pickle.dump(D,open(outfn,'wb'))