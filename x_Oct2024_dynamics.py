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
import xarray as xr

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)


home = '/data2/pmr4/eab32/'

#set dt list based on samples
Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

lag=1 #testing delayed release
first_sample_utc = Pacific.localize(datetime(2024,10,16,10,00)).astimezone(utc)+timedelta(hours=lag)
last_sample_utc = Pacific.localize(datetime(2024,10,16,18,00)).astimezone(utc)+timedelta(hours=lag)
dt_list0 = pd.date_range(start=first_sample_utc,end = last_sample_utc, freq="15min").to_pydatetime().tolist()
ts_list = [datetime.timestamp(dt) for dt in dt_list0]
nt = len(ts_list)

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

subscript = '1945rel_dye'

his_dir0 = '/data1/jxiong/LO_roms/hc12_v00_vldye/'
# from 11am to 6pm on 10/16/2024 PDT
his_dir_list = [his_dir0 + 'f2024.10.16.'+subscript,
                his_dir0 + 'f2024.10.17.'+subscript
            ]
f_list = []
flag=0
print('preprocessing...')
for his_dir in his_dir_list:
    his_dir_fn_list = os.listdir(his_dir)
    fnh = [x for x in his_dir_fn_list if x[:9]=='ocean_his']
    fnh.sort()
    for fn in fnh:
        ds = nc.Dataset(his_dir+'/'+fn)
        ot = ds['ocean_time'][:]
        dt = utc.localize(datetime(1970,1,1)+timedelta(seconds=ot[0]))
        if (dt>first_sample_utc) & (dt<last_sample_utc):
            f_list.append(his_dir+'/'+fn)
        ds.close()

f_list.sort()
# ntt = len(f_list)

ds = xr.open_mfdataset(f_list, combine = 'by_coord')#, concat_dim = 'time')
outfn = '/data2/pmr4/eab32/LO_output/extract/hc12_v00_vldye/ocean_his_HC2_VLdye_DPmask_LOcorr_1hrlag_combined.nc'
ds.to_netcdf(outfn) # Export netcdf file
print('saved to {}'.format(outfn))

