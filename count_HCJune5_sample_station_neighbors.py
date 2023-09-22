import os
import sys
import numpy as np
from datetime import datetime, timedelta
import pytz
import netCDF4 as nc
import pandas as pd
import efun
import pickle

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)

home = '/data2/pmr4/eab32/'

#find June sampling particle tracks
track_dir0 = home+'LO_output/tracks/'

f_list = os.listdir(track_dir0)
f_list.sort()
f_list = [x for x in f_list if (x[:28]=='hc_bangor_perimeter_point_3d')&(x[-4]=='6')]

#set reference time
Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

first_sample_utc = Pacific.localize(datetime(2023,6,5,10,30)).astimezone(utc)
first_release = first_sample_utc-timedelta(days=2)
last_sample_utc = Pacific.localize(datetime(2023,6,5,15,30)).astimezone(utc)
dt_list0 = pd.date_range(start=first_sample_utc,end = last_sample_utc, freq="15min").to_pydatetime().tolist()

#import sampling locations
sampling_data_fn = home+'LO_data/eDNA/June5_sampling_locations.csv'
dfs = pd.read_csv(sampling_data_fn,sep=',',engine='python')
dfs['approx_dt_pdt']=[Pacific.localize(pd.to_datetime(dfs.date_collected.values[i]+' '+dfs.approx_time_collected.values[i],utc=False)) for i in range(len(dfs.date_collected.values))]
dfs['mid_dt_pdt']=[Pacific.localize(pd.to_datetime(dfs.date_collected.values[i]+' '+dfs.mid_time_collected.values[i],utc=False)) for i in range(len(dfs.date_collected.values))]
dfs['approx_dt']=[dt.astimezone(utc) for dt in dfs.approx_dt_pdt]
dfs['mid_dt']=[dt.astimezone(utc) for dt in dfs.mid_dt_pdt]

figcount =0
radius_list = [10,100,1000]
depth_list = [1,2,3,4,5]

nradius = len(radius_list)
ndepth = len(depth_list)
nstation = len(dfs.sampling_location)

station_neighbors = np.zeros((nstation,ndepth,nradius))

# unlike plotting code, open netCDF files only once and loop through timesteps within files
for f in f_list:

    track_dir = track_dir0+f

    # build a keyname from the release filename
    file_list = os.listdir(track_dir)
    file_list = [x for x in file_list if x[:3]=='rel']
    rel_fn = file_list[0]

    ds = nc.Dataset(track_dir+'/'+rel_fn)

    ot = ds['ot'][:].data
    
    dt_list = [pd.Timestamp(utc.localize(datetime(1970,1,1)+timedelta(seconds=ott))) for ott in ot]

    ot0 = dt_list[0]
    ot1 = dt_list[-1]
    
    # loop through stations
    for station in range(nstation):
        
        dt = dfs.approx_dt[station]
        
        # check if sample was taken within the lifetime of the particle release
        if (dt>ot0)&(dt<ot1):
            
            t = argnearest(dt_list,dt)
            
            # loop through depth thresholds
            for depth in range(ndepth):
                
                refdepth = depth_list[depth]
                zmask = ds['z'][t,:]>(ds['zeta'][t,:]-refdepth)
                
                if np.sum(zmask)>0: # don't bother looking at horizontal threshold if no particles are within depth threshold
                    
                    reflon = dfs.long[station]
                    reflat = dfs.lat[station]
                    [x,y] = efun.ll2xy(ds['lon'][t,:],ds['lat'][t,:],reflon,reflat)
                    dist = np.sqrt(x**2+y**2)
                    
                    # loop through horizontal distance thresholds
                    for radius in range(nradius):
                        
                        refradius = radius_list[radius]
                        distmask = dist<refradius
                        
                        # combine horizontal and vertical masks
                        allmask = distmask*zmask
                        allmaskcount = np.sum(allmask)
                        
                        station_neighbors[station,depth,radius]+=allmaskcount
                        
                
    ds.close()
    
# now put in a dict using the field ID of the sampling location
readme = 'station_neighbors indexes [ndepth,nradius]\nBoth radius_list and depth_list are in meters'
D = {}
for station in range(nstation):
    fieldID = dfs.sampling_location[station]
    D[fieldID] = {}
    # first save major result
    D[fieldID]['station_neighbors'] = station_neighbors[station,:]
    
    # now include some metadata - will be same for each fieldID
    for var in ['depth_list','radius_list','readme']:
        D[fieldID][var] = locals()[var]
        
    # include location and time of sample station by fieldID
    for var in ['long','lat','mid_dt','mid_dt_pdt','approx_dt','approx_dt_pdt','eDNA_present']:
        D[fieldID][var] = dfs[var][station]

# finally, pickle
outfn = home+'LO_data/eDNA/HCJune5_bangor_perimeter_point_sampling_station_neighbors.p'
pickle.dump(D,open(outfn,'wb'))
print(f'Saved to {outfn}')
