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
f_list = [x for x in f_list if (x[:11]=='hc_dolph_3d')&(x[17:24]=='2024.06')]
# note current these are the only releases I've done with the new tracker code, so no need to subselect

#set dt list based on samples
Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

first_sample_utc = Pacific.localize(datetime(2024,6,13,10,0)).astimezone(utc)
last_sample_utc = Pacific.localize(datetime(2024,6,13,17,0)).astimezone(utc)
dt_list0 = pd.date_range(start=first_sample_utc,end = last_sample_utc, freq="60min").to_pydatetime().tolist()
ts_list = [datetime.timestamp(dt) for dt in dt_list0]
nt = len(ts_list)

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

#set up heatmap bins
# start by using same as model grid
# grid_fn = track_dir0+f_list[0]+'/grid.nc'
his_dir = '/data1/jxiong/LO_roms/hc11_v01_uu0k/'
his_fn_list = [his_dir + 'f2024.06.13/ocean_his_0018.nc',
                his_dir + 'f2024.06.13/ocean_his_0019.nc',
                his_dir + 'f2024.06.13/ocean_his_0020.nc',
                his_dir + 'f2024.06.13/ocean_his_0021.nc',
                his_dir + 'f2024.06.13/ocean_his_0022.nc',
                his_dir + 'f2024.06.13/ocean_his_0023.nc',
                his_dir + 'f2024.06.13/ocean_his_0024.nc',
                his_dir + 'f2024.06.13/ocean_his_0025.nc']
his_fn_list.sort()
# dsg = nc.Dataset(his_fn_list[0])
# # G,S,T = zrfun.get_basic_info(grid_fn)
# #use psi grid instead of rho to define box edges
# x_edges,y_edges = efun.ll2xy(dsg['lon_psi'][0,:],dsg['lat_psi'][:,0],lon0,lat0)
# nx = len(x_edges)-1
# ny = len(y_edges)-1
# #ugh I can't actually get the z's without zeta, which requires opening the big model output files.
# #z_w = zrfun.get_z(G['h'],)
# dsg.close()

# for now, use regularly shaped zvec. Deepest sample = 121 m
# z_edges = np.linspace(-150,0,16) # produces a vector [-150, -140, -130, ... -20, -10, 0]
# nz = np.shape(dsg['salt'][0,:])[0]-1
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
        S = zrfun.get_basic_info(fn,only_S=True)
    zeta = ds['zeta'][0,:]
    z_w0 = zrfun.get_z(h, zeta, S, only_w=True)
    z_w_x = 0.5*(z_w0[:,:,1:]+z_w0[:,:,:-1])
    z_w_xy = 0.5*(z_w_x[:,1:,:]+z_w_x[:,:-1,:])
    z_edges[tt,:] = z_w_xy[1:-1,:,:]
    ds.close()
    tt+=1

particle_map = np.zeros((nt,nz-2,ny-2,nx-2))
# particle_age_lists = [[[[[] for x in range(nx-2)] for y in range(ny-2)] for z in range(nz-2)] for t in range(nt)]

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
            
            xymask = (yp>y_edges[yi-1])&(yp<y_edges[yi])&(xp>x_edges[xi-1])&(xp<x_edges[xi])
            
            hist,edges = np.histogram(ds['z'][pt,xymask],z_edges[t,:,yi,xi])
            particle_map[t,:,yi,xi] += hist
                
        # hist, edges = np.histogramdd(np.array([xp,yp,ds['z'][pt,:]]),np.array([x_edges,y_edges,z_edges[t,:,:,:]]))
        # particle_map[t,:,:,:] = hist
        # for z in range(nz):
        #     z0 = z_edges[z]
        #     z1 = z_edges[z+1]
        #     zmask = (ds['z'][pt,:]>z0)&(ds['z'][pt,:]<z1)
        #
        #     if np.sum(zmask)>0:
        #         xp,yp = efun.ll2xy(ds['lon'][pt,zmask],ds['lat'][pt,zmask],lon0,lat0)
        #         hist = np.histogram2d(xp,yp,bins=[x_edges,y_edges])[0].T
        #
        #         particle_map[t,z,:]+= hist
                
                # # make lists of particle ages
                # if np.any(hist):
                #     yinds = np.argwhere(np.any(hist,axis=1))[0]
                #     for yi in yinds:
                #         xinds = np.argwhere(hist[yi,:]>0)[0]
                #         if len(xinds)>0:
                #             for xi in xinds:
                #                 age_list = [delta_T]*int(hist[yi,xi])
                #                 for age in age_list:
                #                     particle_age_lists[t][z][yi][xi].append(age)

    ds.close()
    count+=1


D = {}
var_list = ['x_edges','y_edges','z_edges','ts_list','particle_map']
for var in var_list:
    D[var] = locals()[var]


outfn = home+'LO_data/eDNA/June2024_3dhist_zw_no_ages.p'
pickle.dump(D,open(outfn, 'wb'))
print('saved to {}'.format(outfn))
