#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate a filtered time series of a lat-lon moor of ROMS output
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


dir0 = '/data1/jxiong/LO_roms/hc11_v01_uu0kdye/f2023.06.05_test_case' #inside folder are ocean_his_0001.nc - ocean_his_0024.nc

fn_list = os.listdir(dir0)
fn_list.sort()
fn_list = [x for x in fn_list if x[:9]=='ocean_his']

nt = len(fn_list)

dt_list = []
tt = 0
for t in range(nt):

    
    oceanfnh = str(dir0) + '/' + fn_list[t]

    dsh = nc4.Dataset(oceanfnh)

    if tt==0:
        
        #load in lat & lon on u & v grids
        lonr = dsh['lon_rho'][:]
        latr = dsh['lat_rho'][:]
        maskr = dsh['mask_rho'][:]
        
        nyr,nxr = np.shape(lonr)
        
        #initialize array
        dye = np.zeros((nt,nyr,nxr))
        
    
    ot = dsh['ocean_time'][:].data[0]
    dt_list.append(datetime(1970,1,1)+timedelta(seconds=ot))
    
    dye[tt,:] = dsh['dye_01'][0,-1, :]
    
    dsh.close()
    tt = tt+1
    


D = dict()
var_list = ['dt_list','lonr','latr','maskr','dye']

for var in var_list:
    D[var] = locals()[var]


out_dir = '/data2/pmr4/eab32/LO_data/eDNA/'
out_fn = out_dir+'HC_surfdye_0605.p'
pickle.dump(D, open(out_fn, 'wb'))
print('saving to %s' % out_fn)
