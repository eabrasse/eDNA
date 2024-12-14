# -*- coding: utf-8 -*-
"""
modified from Parker's code
(c) Elizabeth Brasseale 11/21/2022

"""
import os
import sys
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
import pytz

utc = pytz.utc

# name output file
dir0='/Users/elizabethbrasseale/Projects/eDNA/data/'

example_fn = dir0+'f2024.10.16_problem/ocean_his_0018.nc'
dsex = nc.Dataset(example_fn)
lonr = dsex['lon_rho'][:]

fn = dir0+'f2024.10.16/LO_vars_for_opendrift.nc'
# get rid of the old version, if it exists
try:
    os.remove(fn)
except OSError:
    pass # assume error was because the file did not exist
ds = nc.Dataset(fn, 'w', format='NETCDF4')

for dname, the_dim in dsex.dimensions.items():
    if 'time' in dname:
        ds.createDimension(dname, 1)
    else:
        ds.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

# add these to new ncfile
vv = ds.createVariable('spherical', int, [])
vv[:] = np.ma.masked_array(data=1,
             mask=False,
       fill_value=999999,
            dtype=int)
vv.units = 'None'
vv.standard_name = 'spherical'      

vv = ds.createVariable('Vtransform', int, [])
vv[:] = np.ma.masked_array(data=2,
             mask=False,
       fill_value=999999,
            dtype=int)
vv.units = 'None'
vv.standard_name = 'Vtransform'  

vv = ds.createVariable('Vstretching', int, [])
vv[:] = np.ma.masked_array(data=4,
             mask=False,
       fill_value=999999,
            dtype=int)
vv.units = 'None'
vv.standard_name = 'Vstretching'     

# GRID ANGLE
vv = ds.createVariable('angle', float, ('eta_rho', 'xi_rho'))
vv[:] = np.ma.masked_array(data=np.zeros(lonr.shape),
             mask=False,
       fill_value=1e+20,
            dtype=float)
vv.units = 'None'
vv.standard_name = 'angle'     

# GLS parameters
vv = ds.createVariable('gls_p', float, [])
vv[:] = np.ma.masked_array(data=3,
             mask=False,
       fill_value=1e+20,
            dtype=float)
vv.units = 'None'
vv.long_name = 'stability exponent'
vv.standard_name = 'gls_p'  

vv = ds.createVariable('gls_m', float, [])
vv[:] = np.ma.masked_array(data=1.5,
             mask=False,
       fill_value=1e+20,
            dtype=float)
vv.units = 'None'
vv.long_name = 'turbulent kinetic energy exponent'
vv.standard_name = 'gls_m' 

vv = ds.createVariable('gls_n', float, [])
vv[:] = np.ma.masked_array(data=-1,
             mask=False,
       fill_value=1e+20,
            dtype=float)
vv.units = 'None'
vv.long_name = 'turbulent length scale exponent'
vv.standard_name = 'gls_n'  

vv = ds.createVariable('gls_Kmin', float, [])
vv[:] = np.ma.masked_array(data=7.6e-6,
             mask=False,
       fill_value=1e+20,
            dtype=float)
vv.units = 'None'
vv.long_name = 'minimum value of specific turbulent kinetic energy'
vv.standard_name = 'gls_Kmin'  

vv = ds.createVariable('gls_Pmin', float, [])
vv[:] = np.ma.masked_array(data=1.0e-12,
             mask=False,
       fill_value=1e+20,
            dtype=float)
vv.units = 'None'
vv.long_name = 'minimum Value of dissipation'
vv.standard_name = 'gls_Pmin' 

vv = ds.createVariable('gls_cmu0', float, [])
vv[:] = np.ma.masked_array(data=0.5477,
             mask=False,
       fill_value=1e+20,
            dtype=float)
vv.units = 'None'
vv.long_name = 'stability coefficient'
vv.standard_name = 'gls_cmu0'   

vv = ds.createVariable('gls_c1', float, [])
vv[:] = np.ma.masked_array(data=1.44,
             mask=False,
       fill_value=1e+20,
            dtype=float)
vv.units = 'None'
vv.long_name = 'shear production coefficient'
vv.standard_name = 'gls_c1'   

vv = ds.createVariable('gls_c2', float, [])
vv[:] = np.ma.masked_array(data=1.92,
             mask=False,
       fill_value=1e+20,
            dtype=float)
vv.units = 'None'
vv.long_name = 'dissipation coefficient'
vv.standard_name = 'gls_c2'  

vv = ds.createVariable('gls_c3m', float, [])
vv[:] = np.ma.masked_array(data=-0.4,
             mask=False,
       fill_value=1e+20,
            dtype=float)
vv.units = 'None'
vv.long_name = 'buoyancy production coefficient (minus)'
vv.standard_name = 'gls_c3m'   

vv = ds.createVariable('gls_c3p', float, [])
vv[:] = np.ma.masked_array(data=1.0,
             mask=False,
       fill_value=1e+20,
            dtype=float)
vv.units = 'None'
vv.long_name = 'buoyancy production coefficient (plus)'
vv.standard_name = 'gls_c3p'

vv = ds.createVariable('gls_sigk', float, [])
vv[:] = np.ma.masked_array(data=1.0,
             mask=False,
       fill_value=1e+20,
            dtype=float)
vv.units = 'None'
vv.long_name = 'constant Schmidt number for TKE'
vv.standard_name = 'gls_sigk'   

vv = ds.createVariable('gls_sigp', float, [])
vv[:] = np.ma.masked_array(data=1.3,
             mask=False,
       fill_value=1e+20,
            dtype=float)
vv.units = 'None'
vv.long_name = 'constant Schmidt number for PSI'
vv.standard_name = 'gls_sigp'   

ds.close()

