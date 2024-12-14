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
# f_list = [x for x in f_list if (x[:11]=='hc_dolph_3d')&(x[-4]!='6')] 
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
ex_his_fn = '/data1/jxiong/LO_roms/hc11_v01_uu0k/f2024.06.13/ocean_his_0002.nc'
dsg = nc.Dataset(ex_his_fn)
# G,S,T = zrfun.get_basic_info(grid_fn)
#use psi grid instead of rho to define box edges
x_edges,y_edges = efun.ll2xy(dsg['lon_psi'][0,:],dsg['lat_psi'][:,0],lon0,lat0)
nx = len(x_edges)-1
ny = len(y_edges)-1
#ugh I can't actually get the z's without zeta, which requires opening the big model output files.
#z_w = zrfun.get_z(G['h'],)

# for now, use regularly shaped zvec. Deepest sample = 121 m
z_edges = np.linspace(-150,0,16) # produces a vector [-150, -140, -130, ... -20, -10, 0]
nz = np.shape(z_edges)[0]-1

dsg.close()

particle_map = np.zeros((nt,nz,ny,nx))
particle_age_lists = [[[[[] for x in range(nx)] for y in range(ny)] for z in range(nz)] for t in range(nt)]

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
        
        for z in range(nz):
            z0 = z_edges[z]
            z1 = z_edges[z+1]
            zmask = (ds['z'][pt,:]>z0)&(ds['z'][pt,:]<z1)
            
            if np.sum(zmask)>0:
                xp,yp = efun.ll2xy(ds['lon'][pt,zmask],ds['lat'][pt,zmask],lon0,lat0)
                hist = np.histogram2d(xp,yp,bins=[x_edges,y_edges])[0].T
                
                particle_map[t,z,:]+= hist
                
                # make lists of particle ages
                if np.any(hist):
                    yinds = np.argwhere(np.any(hist,axis=1))[0]
                    for yi in yinds:
                        xinds = np.argwhere(hist[yi,:]>0)[0]
                        if len(xinds)>0:
                            for xi in xinds:
                                age_list = [delta_T]*int(hist[yi,xi])
                                for age in age_list:
                                    particle_age_lists[t][z][yi][xi].append(age)

    ds.close()
    count+=1


D = {}
var_list = ['x_edges','y_edges','z_edges','ts_list','particle_map','particle_age_lists']
for var in var_list:
    D[var] = locals()[var]


outfn = home+'LO_data/eDNA/June2024_3dhist_w_ages.p'
pickle.dump(D,open(outfn, 'wb'))
print('saved to {}'.format(outfn))
