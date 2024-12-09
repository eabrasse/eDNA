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

hc11 = {'grid':'unmasked'}
hc12 = {'grid':'masked'}
hc11['dir0'] = '/data1/jxiong/LO_roms/hc11_v01_uu0k/'
hc12['dir0'] = '/data1/jxiong/LO_roms/hc12_v00_uu0k/'
for hc in hc11, hc12:
    hc['dt_list'] = []
    hc['f_list'] = [hc['dir0'] + 'f2024.10.16',
        hc['dir0'] + 'f2024.10.17'
                ]
                
    hc['f_list'].sort()

nday = len(hc11['f_list'])
nhr = nday*24+1

vname_list = ['u','v','ubar','vbar']


for hc in hc11, hc12:
    print(hc['grid'])
    tt=0
    fcount = 1
    for f in hc['f_list']:
    
        print('working on folder {} of {}'.format(fcount,len(hc['f_list'])))

        # build a keyname from the release filename
        his_list = os.listdir(f)
        his_list = [x for x in his_list if x[:9]=='ocean_his']
        his_list.sort()
    
        for his_fn in his_list:

            filefn = f+'/'+his_fn
            ds = nc.Dataset(filefn)

            ot = ds['ocean_time'][:].data
            hc['dt_list'].append(datetime(1970,1,1)+timedelta(seconds=ot[0]))
        
            if tt==0:
                hc['lonr'] = ds['lon_rho'][:]
                hc['latr'] = ds['lat_rho'][:]
                hc['maskr'] = ds['mask_rho'][:]
                hc['h'] = ds['h'][:]
            
                for var in vname_list:
                    # build up the shape of the variable
                    # based on the total time steps + grid shape
                    varshape0 = ds[var].shape[:]
                    varshape = [nhr]
                    for nxy in varshape0[1:]:
                        varshape.append(nxy)

                    hc[var] = np.zeros(varshape)
                
            for var in vname_list:
                hc[var][tt,:] = ds[var][0,:]
            
            ds.close()
            tt+=1
        fcount+=1

D = {'hc11':hc11,'hc12':hc12}

outfn = home+'LO_data/eDNA/Oct2024_hc11_hc12_vel_corr.p'

pickle.dump(D,open(outfn, 'wb'))
print('saved to {}'.format(outfn))
