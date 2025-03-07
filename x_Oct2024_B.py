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
#find June sampling particle tracks
track_dir0 = home+'LO_output/tracks2/hc11_v01_uu0k/'

f_list = os.listdir(track_dir0)
f_list.sort()
f_list = [x for x in f_list if (x[:12]=='hc_dolph2_3d')&(x[18:25]=='2024.10')]
# note current these are the only releases I've done with the new tracker code, so no need to subselect

# since particles have consistent release points, we can derive the VL pen only releases
VL_only = True

track_dir = track_dir0+f_list[0]
file_list = os.listdir(track_dir)
file_list = [x for x in file_list if x[:3]=='rel']
rel_fn = file_list[0]

filefn = track_dir+'/'+rel_fn

ds = nc.Dataset(filefn)
if VL_only:
    pmask = ds['lat'][0,:]<47.741
else:
    pmask = np.ones(ds['lat'][0,:].shape)
        
ds.close()

#set dt list based on samples
Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

first_sample_utc = Pacific.localize(datetime(2024,10,16,11,00)).astimezone(utc)
last_sample_utc = Pacific.localize(datetime(2024,10,16,18,00)).astimezone(utc)
dt_list0 = pd.date_range(start=first_sample_utc,end = last_sample_utc, freq="15min").to_pydatetime().tolist()
ts_list = [datetime.timestamp(dt) for dt in dt_list0]
nt = len(ts_list)

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

#set up heatmap bins
# start by using same as model grid
# grid_fn = track_dir0+f_list[0]+'/grid.nc'
his_dir = '/data1/jxiong/LO_roms/hc11_v01_uu0k/'
# from 11am to 6pm on 10/16/2024 PDT
his_fn_list = [his_dir + 'f2024.10.16/ocean_his_0019.nc',
                his_dir + 'f2024.10.16/ocean_his_0020.nc',
                his_dir + 'f2024.10.16/ocean_his_0021.nc',
                his_dir + 'f2024.10.16/ocean_his_0022.nc',
                his_dir + 'f2024.10.16/ocean_his_0023.nc',
                his_dir + 'f2024.10.16/ocean_his_0024.nc',
                his_dir + 'f2024.10.16/ocean_his_0025.nc',
                his_dir + 'f2024.10.17/ocean_his_0002.nc'
            ]
his_fn_list.sort()

ts_list_m = []
tt=0
for fn in his_fn_list:
    ds = nc.Dataset(fn)
    if fn==his_fn_list[0]:
        nz,ny,nx = ds['salt'][0,:,:,:].shape
        x_edges,y_edges = efun.ll2xy(ds['lon_psi'][0,:],ds['lat_psi'][:,0],lon0,lat0)
        # x_edges = np.tile(np.reshape(x_edges,(1,1,nx-1)),(nz-1,ny-1,1))
        # y_edges = np.tile(np.reshape(y_edges,(1,ny-1,1)),(nz-1,1,nx-1))
        z_edges = np.zeros((len(his_fn_list),nz-1,ny-1,nx-1))
        h = ds['h'][:]
        G,S,T = zrfun.get_basic_info(fn)
        lonr = G['lon_rho'][:]
        latr = G['lat_rho'][:]
        xr,yr = efun.ll2xy(lonr,latr,lon0,lat0)
    zeta = ds['zeta'][0,:]
    z_w0 = zrfun.get_z(h, zeta, S, only_w=True)
    z_w_x = 0.5*(z_w0[:,:,1:]+z_w0[:,:,:-1])
    z_w_xy = 0.5*(z_w_x[:,1:,:]+z_w_x[:,:-1,:])
    z_edges[tt,:] = z_w_xy[1:-1,:,:]
    ot = ds['ocean_time'][:]
    ts_list_m.append(datetime.timestamp(datetime(1970,1,1,tzinfo=pytz.utc)+timedelta(seconds=ot[0])))
    ds.close()
    tt+=1
ts_list_m = np.array(ts_list_m)
# #move z_edges time steps from hourly to between hours
# z_edges = 0.5*(z_edges[1:,:]+z_edges[:-1,:])
z_edges_15min = np.zeros((nt,nz-1,ny-1,nx-1))
for tt in range(nt):
    t0m = np.argwhere(ts_list_m<=ts_list[tt])[-1][0]
    t1m = np.argwhere(ts_list_m>=ts_list[tt])[0][0]
    # if we're on a time step, t0m==t1m and z_edges_15min = z_edges[t0m,:] + 0
    # welllll maybe to avoid dividing by zero...
    if t0m==t1m:
        z_edges_15min[tt,:] = z_edges[t0m,:]
    else:
        z_edges_15min[tt,:] = z_edges[t0m,:] + (z_edges[t1m,:]-z_edges[t0m,:])*(ts_list[tt]-ts_list_m[t0m])/(ts_list_m[t1m]-ts_list_m[t0m])
    

particle_map = np.zeros((nt,nz-2,ny-2,nx-2))
# particle_age_lists = [[[[[] for x in range(nx-2)] for y in range(ny-2)] for z in range(nz-2)] for t in range(nt)]

# VLlon = -122.733598
# VLlat = 47.740000
# VLx, VLy = efun.ll2xy(VLlon,VLlat,lon0,lat0)
# VL_particle_profile = np.zeros((nt,nz-2))
# VLxi = np.argmin(np.abs(xr[0,1:-1]-VLx))
# VLyi = np.argmin(np.abs(yr[1:-1,0]-VLy))
#
# HAlat = 47.742236
# HAlon = -122.729975
# HAx, HAy = efun.ll2xy(HAlon,HAlat,lon0,lat0)
# HA_particle_profile = np.zeros((nt,nz-2))
# HAxi = np.argmin(np.abs(xr[0,1:-1]-HAx))
# HAyi = np.argmin(np.abs(yr[1:-1,0]-HAy))

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
VL['profile'] = np.zeros((nt,nz-2))

HA = {}
HA['name'] = 'Husbandry\nArea'
HA['lat'] = 47.742236
HA['lon'] = -122.729975
HA['x'], HA['y'] = efun.ll2xy(HA['lon'],HA['lat'],lon0,lat0)
# HA_particle_profile = np.zeros((nt,nz-2))
HA['xi'] = np.argmin(np.abs(xr[0,1:-1]-HA['x']))
HA['yi'] = np.argmin(np.abs(yr[1:-1,0]-HA['y']))
# HA['col'] = tab10(6)
HA['profile'] = np.zeros((nt,nz-2))

NB = {}
NB['name'] = 'NOAA\nBoat'
NB['lat'] = 47.736613
NB['lon'] = -122.743109
NB['x'],NB['y'] = efun.ll2xy(NB['lon'],NB['lat'],lon0,lat0)
NB['xi'] = np.argmin(np.abs(xr[0,1:-1]-NB['x']))
NB['yi'] = np.argmin(np.abs(yr[1:-1,0]-NB['y']))
# NB['col'] = tab10(9)
NB['profile'] = np.zeros((nt,nz-2))

station_list = [VL,HA,NB]

count=1

for f in f_list:
    
    print(f'working on file {count} of {len(f_list)}')

    track_dir = track_dir0+f

    # build a keyname from the release filename
    file_list = os.listdir(track_dir)
    file_list = [x for x in file_list if x[:3]=='rel']
    rel_fn = file_list[0]

    filefn = track_dir+'/'+rel_fn
    ds = nc.Dataset(filefn)

    ot = ds['ot'][:].data
    ts_list_p = []
    for tt in ot:
        ts_list_p.append(datetime.timestamp(datetime(1970,1,1,tzinfo=pytz.utc)+timedelta(seconds=tt)))
    
    for t in range(nt):
        print(f'Time step {t}...')
        dt_list = [np.abs(ts_p-ts_list[t]) for ts_p in ts_list_p]
        pt = np.argmin(dt_list)
        # pt = np.argmin(np.abs(ts_list_p-ts_list[t]))
        delta_T = ts_list_p[pt]-ts_list_p[0]
        
        xp,yp = efun.ll2xy(ds['lon'][pt,:],ds['lat'][pt,:],lon0,lat0)
        
        xis = np.digitize(xp,x_edges)
        yis = np.digitize(yp,y_edges)
        
        xyi = [[xis[i],yis[i]] for i in range(len(xis))]
        
        xyi_u = np.unique(xyi,axis=0)
        
        for [xi,yi] in xyi_u:
            
            if (xi in [0,len(x_edges)]) or (yi in [0,len(y_edges)]):
                continue
            
            xymask = (yp>y_edges[yi-1])&(yp<y_edges[yi])&(xp>x_edges[xi-1])&(xp<x_edges[xi])&pmask
            
            hist,edges = np.histogram(ds['z'][pt,xymask],z_edges_15min[t,:,yi,xi])
            particle_map[t,:,yi-1,xi-1] += hist
            
        
        for station in station_list:
            station['rpm'] = (np.sqrt((xp-station['x'])**2+(yp-station['y'])**2)<100)&pmask
            hist,edges = np.histogram(ds['z'][pt,station['rpm']],z_edges_15min[t,:,station['yi'],station['xi']])
            station['profile'][t,:] += hist


    ds.close()
    count+=1


D = {}
var_list = ['x_edges','y_edges','z_edges_15min','ts_list','particle_map','station_list']
for var in var_list:
    D[var] = locals()[var]


outfn = home+'LO_data/eDNA/Oct2024_3dhist_zw_no_ages_15min'
if VL_only:
    outfn = outfn + "_VL_only.p"
else:
    outfn = outfn + '.p'

pickle.dump(D,open(outfn, 'wb'))
print('saved to {}'.format(outfn))
