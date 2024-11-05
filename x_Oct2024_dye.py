import os
import sys
import numpy as np
from datetime import datetime, timedelta
import pytz
import netCDF4 as nc
import pandas as pd
import pickle
import efun
from lo_tools import zrfun

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)


home = '/data2/pmr4/eab32/'

#set dt list based on samples
Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

first_sample_utc = Pacific.localize(datetime(2024,10,16,10,00)).astimezone(utc)
last_sample_utc = Pacific.localize(datetime(2024,10,16,18,00)).astimezone(utc)
dt_list0 = pd.date_range(start=first_sample_utc,end = last_sample_utc, freq="15min").to_pydatetime().tolist()
ts_list = [datetime.timestamp(dt) for dt in dt_list0]
nt = len(ts_list)

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

subscript = '1845rel'
#set up heatmap bins
# start by using same as model grid
# grid_fn = track_dir0+f_list[0]+'/grid.nc'
his_dir0 = '/data1/jxiong/LO_roms/hc11_v01_vldye/'
# from 11am to 6pm on 10/16/2024 PDT
his_dir_list = [his_dir0 + 'f2024.10.16.'+subscript,
                his_dir0 + 'f2024.10.17.'+subscript
            ]
f_list = []
flag=0
print('preprocessing...')
for his_dir in his_dir_list:
    his_dir_fn_list = os.listdir(his_dir)
    fnh = [x for x in his_dir_fn_list if x[:9]=='ocean_his']
    fnh.sort()
    for fn in fnh:
        ds = nc.Dataset(his_dir+'/'+fn)
        if flag==0:
            nz,ny,nx = ds['salt'][0,:,:,:].shape
            lonr = ds['lon_rho'][:]
            latr = ds['lat_rho'][:]
            xr,yr = efun.ll2xy(lonr,latr,lon0,lat0)
            h = ds['h'][:]
            S= zrfun.get_basic_info(his_dir+'/'+fn, only_S=True)
            flag+=1
        ot = ds['ocean_time'][:]
        dt = utc.localize(datetime(1970,1,1)+timedelta(seconds=ot[0]))
        if (dt>first_sample_utc) & (dt<last_sample_utc):
            f_list.append(his_dir+'/'+fn)
        ds.close()

f_list.sort()
ntt = len(f_list)




# VL pen and Husbandry area sampling stations
VL = {}
VL['name'] = 'VL\nPen\nSouth'
VL['lon'] = -122.733598
VL['lat'] = 47.740000
VL['x'], VL['y'] = efun.ll2xy(VL['lon'],VL['lat'],lon0,lat0)
# VL_particle_profile = np.zeros((nt,nz-2))
VL['xi'] = np.argmin(np.abs(xr[0,1:-1]-VL['x']))
VL['yi'] = np.argmin(np.abs(yr[1:-1,0]-VL['y']))
# VL['col'] = tab10(8)
VL['profile'] = np.zeros((nt,nz))
VL['zr'] = np.zeros((nt,nz))

HA = {}
HA['name'] = 'Husbandry\nArea'
HA['lat'] = 47.742236
HA['lon'] = -122.729975
HA['x'], HA['y'] = efun.ll2xy(HA['lon'],HA['lat'],lon0,lat0)
# HA_particle_profile = np.zeros((nt,nz-2))
HA['xi'] = np.argmin(np.abs(xr[0,1:-1]-HA['x']))
HA['yi'] = np.argmin(np.abs(yr[1:-1,0]-HA['y']))
# HA['col'] = tab10(6)
HA['profile'] = np.zeros((nt,nz))
HA['zr'] = np.zeros((nt,nz))

NB = {}
NB['name'] = 'NOAA\nBoat'
NB['lat'] = 47.736613
NB['lon'] = -122.743109
NB['x'],NB['y'] = efun.ll2xy(NB['lon'],NB['lat'],lon0,lat0)
NB['xi'] = np.argmin(np.abs(xr[0,1:-1]-NB['x']))
NB['yi'] = np.argmin(np.abs(yr[1:-1,0]-NB['y']))
# NB['col'] = tab10(9)
NB['profile'] = np.zeros((nt,nz))
NB['zr'] = np.zeros((nt,nz))

station_list = [VL,HA,NB]
dye_upper10m = np.zeros((nt,ny,nx))
dt_list = []

tt = 0

for f in f_list:
    
    print(f'extracting dye from file {tt} of {len(f_list)}')

    ds = nc.Dataset(f)

    ot = ds['ocean_time'][:]
    dt_list.append((datetime(1970,1,1,tzinfo=pytz.utc)+timedelta(seconds=ot[0])).astimezone(utc))
    
    zeta = ds['zeta'][0,:]
    zr = zrfun.get_z(h, zeta, S, only_r=True)
    
    zrmask = zr[0,:,:,:]<-10 # mask everything deeper than 10m
    dye_masked = np.ma.array(ds['dye_01'][0,:,:,:],zrmask)
    dye_upper10m[tt,:,:] = np.mean(dye_masked,axis=0)
    
    for station in station_list:
        station['zr'][tt,:] = zr[0,:,yi,xi]
        station['profile'][tt,:] = ds['dye_01'][0,:,yi,xi]
    
    ds.close()
    tt+=1



D = {}
var_list = ['xr','yr','dt_list','dye_upper10m','station_list']
for var in var_list:
    D[var] = locals()[var]


outfn = home+'LO_data/eDNA/Oct2024_'+subscript+'_dye_15min.p'


pickle.dump(D,open(outfn, 'wb'))
print('saved to {}'.format(outfn))
