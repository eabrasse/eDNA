#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate a filtered time series of a lat-lon slice of ROMS output
"""

# setup
import os
import sys
from lo_tools import Lfun, zrfun, zfun
Ldir = Lfun.Lstart()

import netCDF4 as nc4
import numpy as np
from datetime import datetime, timedelta
import pickle

gridname = 'hc11'
tag = 'v1'
ex_name = 'uu0k'
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
dir0 = Ldir['roms_out'] / Ldir['gtagex']
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex']
Lfun.make_dir(out_dir)

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[:5]=='f2023']

f_list = f_list[4:9]

ndays = len(f_list)
nt = ndays*24

dt_list = []
tt = 0
for t in range(ndays):
    print('day %i' % t)
    
    fn_list = os.listdir(str(dir0) + '/' + f_list[t])
    fn_list.sort()
    fnh_list = [x for x in fn_list if x[6:9]=='his']
    
    hrs = len(fnh_list)
    for hr in range(hrs):
        
        oceanfnh = str(dir0) + '/' + f_list[t] + '/' + fnh_list[hr]
    
        dsh = nc4.Dataset(oceanfnh)
    
        if tt==0:
            mask_rho = dsh['mask_rho'][:]
            NYr, NXr = np.shape(mask_rho)

            wetdry_mask_rho = np.zeros((nt,NYr,NXr))
            zeta = np.zeros((nt,NYr,NXr))
            
        #dump data slices into filterable matrices
        wetdry_mask_rho[tt,:,:] = dsh['wetdry_mask_rho'][0,:,:]
        zeta[tt,:,:] = dsh['zeta'][0,:,:]
        
        ot = dsh['ocean_time'][:].data[0]
        dt_list.append(datetime(1970,1,1)+timedelta(seconds=ot))
        
        dsh.close()
        tt = tt+1


D = dict()

D['dt_list'] = dt_list
D['zeta'] = zeta
D['mask_rho'] = mask_rho
D['wetdry_mask_rho'] = wetdry_mask_rho

out_fn = out_dir / 'wetdry_mask_tides.p'
pickle.dump(D, open(out_fn, 'wb'))
print('saving to %s' % out_fn)
