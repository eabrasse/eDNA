import os
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from PIL import Image
import numpy as np
from datetime import datetime, timedelta
import pytz
import netCDF4 as nc
import pandas as pd
import pickle

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)

plt.close('all')
home = '/data2/pmr4/eab32/'

#find June sampling particle tracks
track_dir0 = home+'LO_output/tracks/'

f_list = os.listdir(track_dir0)
f_list.sort()
f_list = [x for x in f_list if (x[:11]=='hc_dolph_3d')&(x[-4]!='6')] # because releases were in Jan & Feb

# load eDNA concentrations
eDNA_fn = home+'LO_data/eDNA/ESP_Feb2023_hourly.csv'
df = pd.read_csv(eDNA_fn,sep=',',engine='python')
df.date = pd.to_datetime(df.date)

C0 = {}
count=0
for f in f_list:

    track_dir = track_dir0+f

    # build a keyname from the release filename
    file_list = os.listdir(track_dir)
    file_list = [x for x in file_list if x[:3]=='rel']
    rel_fn = file_list[0]

    ds = nc.Dataset(track_dir+'/'+rel_fn)

    ot = ds['ot'][:].data
    T0 = datetime(1970,1,1,tzinfo=pytz.utc)+timedelta(seconds=ot[0]) #datetime of release
    T1 = datetime(1970,1,1,tzinfo=pytz.utc)+timedelta(seconds=ot[-1])
    if count==0:
        Tmin = T0
        Tmax = T1
    else:
        Tmin = min([Tmin,T0])
        Tmax = max([Tmax,T1])
    C0[rel_fn] = df[df.date==T0].PB_quantity_mean.values[0] 
    
    count+=1

#set reference time
Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

# first_sample_utc = Pacific.localize(datetime(2023,1,31,0,0)).astimezone(utc)
# first_release = first_sample_utc-timedelta(days=2)
# last_sample_utc = Pacific.localize(datetime(2023,2,2,0,0)).astimezone(utc)
# might... do this different
dt_list0 = pd.date_range(start=Tmin,end = Tmax, freq="60min").to_pydatetime().tolist()

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

moor0 = {}
moor0['lon'] = lon0
moor0['lat'] = lat0
moor0['label'] = 'Dolphin pen'

moor1 = {}
moor1['lon'] = -122.72845
moor1['lat'] = 47.747933
moor1['label'] = 'Marginal Pier'

moor2 = {}
moor2['lon'] = -122.729167
moor2['lat'] = 47.745633
moor2['label'] = 'near Marginal Pier'

moor3 = {}
moor3['lon'] = -122.728933
moor3['lat'] = 47.74325
moor3['label'] = 'near Delta Pier'

moor4 = {}
moor4['lon'] = -122.729783
moor4['lat'] = 47.742667
moor4['label'] = 'Delta Pier'

moor_list = [moor0,moor1,moor2,moor3,moor4]


# define domain & bins as done in previous heatmap plotting code
# plot region around delta pier, using release point as indicator
pad = .01
aa = [lon0-pad, lon0+pad,
   lat0-pad, lat0+pad]
zref = -2

nbins = 100
bin_lon_edges=np.linspace(aa[0], aa[1],nbins+1)
bin_lat_edges=np.linspace(aa[0], aa[1],nbins+1)
xx, yy = np.meshgrid(bin_lon_edges[:-1]+0.5*(bin_lon_edges[1]-bin_lon_edges[0]),bin_lat_edges[:-1]+0.5*(bin_lat_edges[1]-bin_lat_edges[0]))

# find which bin
for moor in moor_list:
    #first bin they're greater than
    moor['lon_bin'] = np.argwhere(moor['lon']>bin_lon_edges)[0][0]
    moor['lat_bin'] = np.argwhere(moor['lat']>bin_lat_edges)[0][0]

nt = len(dt_list0)
const_particle_bin = np.zeros((nt,nbins,nbins))
TV_particle_bin = np.zeros((nt,nbins,nbins))

k_decay = 0.02/3600 #units: data 0.02 1/hr, multiply by hr/sec to get 1/sec

count =0
for dt in dt_list0:
    
    #FIGURE OUT HOW TO INCORPORATE TIME VARYING WITHOUT LOOKING IT UP OVER AND OVER
    for f in f_list:

        track_dir = track_dir0+f

        # build a keyname from the release filename
        file_list = os.listdir(track_dir)
        file_list = [x for x in file_list if x[:3]=='rel']
        rel_fn = file_list[0]

        ds = nc.Dataset(track_dir+'/'+rel_fn)

        ot = ds['ot'][:].data

        ot0 = utc.localize(datetime(1970,1,1)+timedelta(seconds=ot[0]))
        ot1 = utc.localize(datetime(1970,1,1)+timedelta(seconds=ot[-1]))

        if (dt>ot0)&(dt<ot1):
            dt_list = [utc.localize(datetime(1970,1,1)+timedelta(seconds=ott)) for ott in ot]
            t = argnearest(dt_list,dt)
            delta_T = ot[t]-ot[0]
            decay = np.exp(-k_decay*delta_T)
            zmask = ds['z'][t,:]>(ds['zeta'][t,:]-2)
            if np.sum(zmask)>0:
                # ax.scatter(ds['lon'][t,zmask],ds['lat'][t,zmask],c='w',s=1,alpha=0.05)
                hist = np.histogram2d(ds['lon'][t,zmask],ds['lat'][t,zmask],bins=[bin_lon_edges,bin_lat_edges])
                
                const_particle_bin[count,:] += decay*hist[0].T
                TV_particle_bin[count,:] += C0[rel_fn]*decay*hist[0].T
                
        ds.close()
    count+=1


for moor in moor_list:
    moor['const_particle_bin'] = constant_particle_bin[:,moor['lat_bin'],moor['lon_bin']]
    moor['TV_particle_bin'] = TV_particle_bin[:,moor['lat_bin'],moor['lon_bin']]
    

outfn = home+'LO_data/eDNA/Feb2023_moorings.p'
pickle.dump(moor,outfn)
