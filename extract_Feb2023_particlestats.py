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
import efun



plt.close('all')
home = '/data2/pmr4/eab32/'

#find June sampling particle tracks
track_dir0 = home+'LO_output/tracks/'

f_list = os.listdir(track_dir0)
f_list.sort()
f_list = [x for x in f_list if (x[:11]=='hc_dolph_3d')&(x[-4]!='6')] # because releases were in Jan & Feb

# dolphin pen location - just need some lon0/lat0 pair to calculate ll2xy
lon0 = -122.729779; lat0 = 47.742773

D = {}

print('Preprocessing particle release experiments')
for f in f_list:
    D[f] = {}

    track_dir = track_dir0+f

    # build a keyname from the release filename
    file_list = os.listdir(track_dir)
    file_list = [x for x in file_list if x[:3]=='rel']
    rel_fn = file_list[0]

    filefn = track_dir+'/'+rel_fn
    D[f]['filefn']=filefn
    ds = nc.Dataset(filefn)
    
    x,y = efun.ll2xy(ds['lon'][:].data,ds['lat'][:].data,lon0,lat0)
    
    # particle cloud COM as a function of time
    D[f]['x_mean'] = np.mean(x,axis=1)
    D[f]['y_mean'] = np.mean(y,axis=1)
    D[f]['z_mean'] = np.mean(ds['z'][:].data,axis=1)
    
    Dx = x-np.tile(np.reshape(D[f]['x_mean'],(D[f]['x_mean'].shape[0],1)),(1,x.shape[1]))
    Dy = y-np.tile(np.reshape(D[f]['y_mean'],(D[f]['y_mean'].shape[0],1)),(1,y.shape[1]))
    Dz = ds['z'][:].data-np.tile(np.reshape(D[f]['z_mean'],(D[f]['z_mean'].shape[0],1)),(1,ds['z'][:].data.shape[1]))
    
    D[f]['avg_dist_from_COM'] = np.mean(np.sqrt(Dx**2+Dy**2+Dz**2),axis=1)

    D[f]['t'] = ds['ot'][:].data-ds['ot'][0].data


outfn = home+'LO_data/eDNA/Feb2023_particlepositionstats.p'
pickle.dump(D,open(outfn, 'wb'))
print('saved to {}'.format(outfn))
