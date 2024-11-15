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

home = '/data2/pmr4/eab32/'

dir0 = '/data1/jxiong/LO_roms/hc11_v01_uu0k/'
f_list = [dir0 + 'f2024.10.14',
    dir0 + 'f2024.10.15',
    dir0 + 'f2024.10.16',
    dir0 + 'f2024.10.17',
    dir0 + 'f2024.10.18'
            ]
f_list.sort()

nday = len(f_list)
nhr = nday*24

# # VL pen and Husbandry area sampling stations
# Brem = {}
# Brem['name'] = 'Bremerton'
# Brem['lat'] = 47.562
# Brem['lon'] = -122.626
#
# Bang = {}
# Bang['name'] = 'Bangor'
# Bang['lat'] = 47.7578
# Bang['lon'] = -122.7306
#
# station_list = [Brem,Bang]

vname_list = ['Uwind','Vwind','svstr','sustr']
vdict = {}
# for station in station_list:
#     for var in var_list:
#         station[var] =  np.zeros((nhr))

dt_list = []

tt=0
fcount = 1
for f in f_list:
    
    print(f'working on folder {fcount} of {len(f_list)}')

    # build a keyname from the release filename
    his_list = os.listdir(f)
    his_list = [x for x in his_list if x[:9]=='ocean_his']
    his_list.sort()
    
    for his_fn in his_list:

        filefn = f+'/'+his_fn
        ds = nc.Dataset(filefn)

        ot = ds['ocean_time'][:].data
        dt_list.append(datetime(1970,1,1)+timedelta(seconds=ot[0]))
        
        if tt==0:
            lonr = ds['lon_rho'][:]
            latr = ds['lat_rho'][:]
            maskr = ds['mask_rho'][:]
            h = ds['h'][:]
            
            for var in vname_list:
                # build up the shape of the variable
                # based on the total time steps + grid shape
                varshape0 = ds[var].shape[:]
                varshape = [nhr]
                for nxy in varshape0[1:]:
                    varshape.append(nxy)
                    
                vdict[var][:] = np.zeros(varshape)
                
        for var in vname_list:
            vdict[var][tt,:] = ds[var][0,:]
            
        ds.close()
        tt+=1
    fcount+=1


othervars = ['dt_list','lonr','latr','h','maskr']
for ovar in othervars:
    vdict[ovar] = locals()[ovar]


outfn = home+'LO_data/eDNA/Oct2024_HC2wind.p'

pickle.dump(vdict,open(outfn, 'wb'))
print('saved to {}'.format(outfn))
