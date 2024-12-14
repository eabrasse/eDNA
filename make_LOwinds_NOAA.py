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
dir0='/Users/elizabethbrasseale/Projects/eDNA/'
# hourly wind file
obs_fn = dir0+'data/IOOS_Wind_Bremerton9445958_20241015_20241018_hourly.csv'
dfo = pd.read_csv(obs_fn,sep=',',engine='python')
dfo['dt'] = pd.to_datetime(dfo['dt_list'],utc=True) #parse datetime string to datetime datatype
# dfo['dt'] = [dt.astimezone(utc) for dt in dfo['dt']] # add timezone awareness
dfo['timestamp'] = dfo['dt'].values.view('int64') // 10 ** 9 # convert from datetime to timestamp

#example file
dir1 = dir0+'data/f2024.10.16_masked/'
# build a keyname from the release filename
his_list = os.listdir(dir1)
his_list = [x for x in his_list if x[:9]=='ocean_his']
his_list.sort()

for his_fn in his_list:

    orig_fn = dir0+'data/f2024.10.16_masked/'+his_fn

    # orig_fn = dir0+'data/f2024.10.16/ocean_his_0016.nc'
    ds1 = nc.Dataset(orig_fn, mode='r')

    #output file
    new_fn = dir0+'data/f2024.10.16_masked_NOAAwinds/'+his_fn

    # get rid of the old version, if it exists
    try:
        os.remove(new_fn)
    except OSError:
        pass # assume error was because the file did not exist
    # ds2 = nc.Dataset(ini_fn, 'w', format='NETCDF3_64BIT_OFFSET')
    ds2 = nc.Dataset(new_fn, 'w', format='NETCDF4')

    # Copy dimensions
    for dname, the_dim in ds1.dimensions.items():
        if 'time' in dname:
            ds2.createDimension(dname, None)
        else:
            ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)


    # copy all variables for forcing file
    for v_name, varin in ds1.variables.items():
        if v_name not in ['Uwind','Vwind']:
            outVar = ds2.createVariable(v_name, varin.datatype, varin.dimensions)
            # Copy variable attributes, {} is a dict comprehension, cool!
            outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
            outVar[:] = varin[:]
            # if varin.ndim > 1:
            #     outVar[:] = varin[0,:]
            # else:
            #     outVar[:] = varin[0]

    # find Uwind and Vwind data that matches timestamp of his file
    ot = ds1['ocean_time'][:]
    #might be the right time stamp already, but just in case... convert to tz-aware datetime and back
    dt = utc.localize(datetime(1970,1,1)+timedelta(seconds=ot[0]))
    ts = dt.timestamp() # convert from datetime to timestamp
    tind = np.argmin(np.abs(dfo['timestamp']-ts))
    print(tind)
    print(dfo.iloc[tind]['dt'])
    print('Uwind = {} ms-1'.format(dfo.iloc[tind].Uwind))
    print('Vwind = {} ms-1'.format(dfo.iloc[tind].Vwind))

    # now make arrays of the necessary shape for the netCDF file,
    # but with a spatially uniform value of the observed wind at Bremerton
    Uwind = np.ones(ds1['Uwind'].shape)*dfo.iloc[tind].Uwind
    Vwind = np.ones(ds1['Vwind'].shape)*dfo.iloc[tind].Vwind

    # add these to new ncfile
    vv = ds2.createVariable('Uwind', float, ds1['Uwind'].dimensions)
    vv[:] = Uwind[:]

    vv = ds2.createVariable('Vwind', float, ds1['Vwind'].dimensions)
    vv[:] = Vwind[:]


    ds1.close()
    ds2.close()

