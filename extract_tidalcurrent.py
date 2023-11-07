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


#locations
tidegauge = {}
tidegauge['name']='TideGauge_9445133'
tidegauge['lon0']=-(122 + 43.8/60)
tidegauge['lat0']=47 + 45.5/60

dolphpen = {}
dolphpen['name']='BangorDolphinPen'
dolphpen['lon0']=-122.729779
dolphpen['lat0']=47.742773

neardolphpen = {}
neardolphpen['name']='NearBangorDolphinPen'
neardolphpen['lon0']=-122.729779
neardolphpen['lat0']=47.742773


# choose which location to extract
print(90*'*')
print('\n%s\n' % '** Choose extract location (return for Tide Gauge) **')
st_list = [tidegauge, dolphpen, neardolphpen]
Nst = len(st_list)
st_dict = dict(zip(range(Nst), st_list))
for nst in range(Nst):
    print(str(nst) + ': ' + st_list[nst]['name'])
my_nst = input('-- Input number -- ')
if len(my_nst)==0:
    location = tidegauge
else:
    location = st_dict[int(my_nst)]

lat_nudge = 0
lon_nudge = 0
if location==dolphpen:
    lat_nudge=-2
if location==neardolphpen:
    lat_nudge=-3
    lon_nudge=3

dt_list = []
tt = 0
for t in range(ndays):
    
    fn_list = os.listdir(str(dir0) + '/' + f_list[t])
    fn_list.sort()
    fnh_list = [x for x in fn_list if x[6:9]=='his']
    
    hrs = len(fnh_list)
    for hr in range(hrs):
        
        oceanfnh = str(dir0) + '/' + f_list[t] + '/' + fnh_list[hr]
    
        dsh = nc4.Dataset(oceanfnh)
    
        if tt==0:

            #initialize arrays
            zeta = np.zeros((nt))
            u = np.zeros((nt))
            v = np.zeros((nt))
            
            #load in lat & lon on u & v grids
            lonr = dsh['lon_rho'][:]
            latr = dsh['lat_rho'][:]
            maskr = dsh['mask_rho'][:]
            lonvecr = lonr[0,:]
            latvecr = latr[:,0]
            
            lonu = dsh['lon_u'][:]
            latu = dsh['lat_u'][:]
            masku = dsh['mask_u'][:]
            lonvecu = lonu[0,:]
            latvecu = latu[:,0]
            
            lonv = dsh['lon_v'][:]
            latv = dsh['lat_v'][:]
            maskv = dsh['mask_v'][:]
            lonvecv = lonv[0,:]
            latvecv = latv[:,0]
            
            # get tidal location index          
            i_r = zfun.find_nearest_ind(lonvecr,location['lon0'])+lon_nudge
            j_r = zfun.find_nearest_ind(latvecr,location['lat0'])+lat_nudge
              
            i_u = zfun.find_nearest_ind(lonvecu,location['lon0'])+lon_nudge
            j_u = zfun.find_nearest_ind(latvecu,location['lat0'])+lat_nudge
            
            i_v = zfun.find_nearest_ind(lonvecv,location['lon0'])+lon_nudge
            j_v = zfun.find_nearest_ind(latvecv,location['lat0'])+lat_nudge #avoid delta pier
    
        
        ot = dsh['ocean_time'][:].data[0]
        dt_list.append(datetime(1970,1,1)+timedelta(seconds=ot))
        
        zeta[tt] = dsh['zeta'][0,j_r,i_r]
        u[tt] = dsh['ubar'][0,j_u,i_u]
        v[tt] = dsh['vbar'][0,j_v,i_v]
        
        
        dsh.close()
        tt = tt+1
        
    print('Finished {:} (day {:} of {:})'.format(dt_list[tt-1].strftime('%b %-d'),(t+1),ndays))


D = dict()
var_list = ['dt_list','zeta','u','v','lonu','latu','masku','lonv','latv','maskv','lonr','latr','maskr']
for var in var_list:
    D[var] = locals()[var]
# D['dt_list'] = dt_list
# D['zeta'] = zeta
# D['u'] = u
# D['v'] = v
D['lon0'] = location['lon0']
D['lat0'] = location['lat0']
D['location_name'] = location['name']


out_fn = out_dir / 'HC_dolph_tidal_current_AT{:}.p'.format(location['name'])
pickle.dump(D, open(out_fn, 'wb'))
print('saving to %s' % out_fn)
