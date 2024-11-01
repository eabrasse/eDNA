#!/usr/bin/env python
"""
Back and forth
==============
"""

import os
import opendrift
# from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
# from opendrift.models.oceandrift import OceanDrift
from opendrift.models.LO_Drift import LO_Drift
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import netCDF4 as nc
import matplotlib
import cmocean as cmo
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import pandas as pd
import pytz

plt.close('all')

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

newcmap = cmo.tools.crop_by_percent(cmo.cm.deep_r, 10, which='max', N=None)

dir0= '/Users/elizabethbrasseale/Projects/eDNA/data/'

meta_fn = dir0+'Drifter metadata.csv'
meta_df = pd.read_csv(meta_fn,sep=',',engine='python')
meta_df['dt0'] = [datetime(2024,10,16,int(meta_df['Release time (PDT)'].values[i][:2]),int(meta_df['Release time (PDT)'].values[i][3:])) for i in range(len(meta_df.index))]
meta_df['dt0'] = meta_df['dt0'].dt.tz_localize('America/Vancouver').dt.tz_convert(pytz.utc)
meta_df['dt1'] = [datetime(2024,10,16,int(meta_df['Recovery time (PDT)'].values[i][:2]),int(meta_df['Recovery time (PDT)'].values[i][3:])) for i in range(len(meta_df.index))]
meta_df['dt1'] = meta_df['dt1'].dt.tz_localize('America/Vancouver').dt.tz_convert(pytz.utc)

relnum = 5
lag = 1

o = LO_Drift(loglevel=20)  # Set loglevel to 0 for debug information


ncfile = dir0+'opendrift/HCdrifter_R{}_{}hrlag.nc'.format(relnum,lag)


try:
    os.remove(ncfile)
except OSError:
    pass
    
ROMS_variable_mapping = {
            # Removing (temoprarily) land_binary_mask from ROMS-variables,
            # as this leads to trouble with linearNDFast interpolation
            'mask_rho': 'land_binary_mask',
            # 'mask_psi': 'land_binary_mask',  # don't want two variables mapping together - raises error now
            'h': 'sea_floor_depth_below_sea_level',
            'zeta': 'sea_surface_height',
            'u': 'x_sea_water_velocity',
            'v': 'y_sea_water_velocity',
            #'u_eastward': 'x_sea_water_velocity',  # these are wrognly rotated below
            #'v_northward': 'y_sea_water_velocity',
            'w': 'upward_sea_water_velocity',
            'temp': 'sea_water_temperature',
            'salt': 'sea_water_salinity',
            'uice': 'sea_ice_x_velocity',
            'vice': 'sea_ice_y_velocity',
            'aice': 'sea_ice_area_fraction',
            'hice': 'sea_ice_thickness',
            'gls': 'turbulent_generic_length_scale',
            'tke': 'turbulent_kinetic_energy',
            'AKs': 'ocean_vertical_diffusivity',
            'sustr': 'surface_downward_x_stress',
            'svstr': 'surface_downward_y_stress',
            'tair': 'air_temperature',
            'wspd': 'wind_speed',
            'uwnd': 'x_wind',
            'vwnd': 'y_wind',
            'uwind': 'x_wind',
            'vwind': 'y_wind',
            'Uwind': 'x_wind',
            'Vwind': 'y_wind',
            }



LO_native = reader_ROMS_native.Reader('/Users/elizabethbrasseale/Projects/eDNA/data/f2024.10.16/ocean_his_*.nc')#,standard_name_mapping=ROMS_variable_mapping)

o.add_reader(LO_native)
# Jilian's options
o.set_config('general:coastline_action', 'previous')  # land boundary
o.set_config('drift:advection_scheme','runge-kutta4')
o.set_config('drift:vertical_mixing',  True)
o.set_config('vertical_mixing:diffusivitymodel', 'environment')
o.set_config('drift:vertical_advection',  True)
o.set_config('general:use_auto_landmask', False)  # use ROMS landmask
o.set_config('general:seafloor_action', 'previous')  # default: lift_to_seafloor (elements are lifted to sea floor level) 
# 
# 
o.disable_vertical_motion()
#%%
# Forward run
# Seeding some particles
#47.742212°N 122.729401°W
lon = meta_df['Release Lon'].loc[relnum]; lat = meta_df['Release Lat'].loc[relnum];
# time = datetime(2024, 10, 16, 19, 10, 0) #18:10 = 11:10am PDT, Husbandry area drifter release
time = meta_df['dt0'].loc[relnum].to_pydatetime().replace(tzinfo=None)+timedelta(hours=int(lag))

# time = LO_native.start_time
nseed = 10000
o.seed_elements(lon, lat, z=0, radius=20, number=nseed, time=time)

td = meta_df['dt1'].loc[relnum]-meta_df['dt0'].loc[relnum]
nsteps = int(np.ceil(td.seconds/150))
# o.run(steps=5*4-1, time_step=-900, time_step_output=900, outfile=ncfile0)
o.run(steps=nsteps, time_step=150, time_step_output=150, outfile=ncfile,
     export_variables=['lon','lat','z','sea_water_temperature','sea_water_salinity','sea_floor_depth_below_sea_level']) #every 2.5 minutes for an hour

#%%
# Print and plot results
print(o)
# o.plot(buffer=.2, fast=True)
# fig = plt.figure(figsize=(8,8))
# ax = fig.gca()
# lons,lats = o.get_lonlats()
# ax.contour(LO_native.lon,LO_native.lat,LO_native.mask_rho.data[0,:,:],colors='gray',linestyles='dashed')
# p=ax.pcolormesh(LO_native.lon,LO_native.lat,LO_native.sea_floor_depth_below_sea_level,cmap=cmo.cm.deep,alpha=0.8)
# for i in range(nseed):
#     ax.plot(lons.data[i,:],lats.data[i,:],color='k',alpha=1,lw=0.3)
# ax.scatter(lons.data[:,0],lats.data[:,0],c='green')
# ax.scatter(lons.data[:,-1],lats.data[:,-1],c='red')
# fudge=0.002
# ax.axis([lons.min()-fudge,lons.max()+fudge,lats.min()-fudge,lats.max()+fudge])
# plt.colorbar(p)
# plt.show(block=False)
# plt.pause(0.1)
#
# #plot particle track depth
# fig = plt.figure(figsize=(8,8))
# ax = fig.gca()
# lons,lats = o.get_lonlats()
#
# norm = plt.Normalize(o.history['z'][:].min(), o.history['z'][:].max())
# for i in range(nseed):
#     points = np.array([lons.data[i,:], lats.data[i,:]]).T.reshape(-1, 1, 2)
#     segments = np.concatenate([points[:-1], points[1:]], axis=1)
#     lc = LineCollection(segments, cmap=newcmap, norm=norm,zorder=10)
#     lc.set_array(o.history['z'][i,:])
#     lc.set_linewidth(0.5)
#     line = ax.add_collection(lc)
# ax.scatter(lons.data[:,0],lats.data[:,0],c='green',zorder=50)
# ax.scatter(lons.data[:,-1],lats.data[:,-1],c=o.history['z'][:,-1],edgecolors='k',zorder=51,cmap=newcmap, norm=norm)
# ax.contour(LO_native.lon,LO_native.lat,LO_native.mask_rho,colors='gray',linestyles='dashed',zorder=200)
# ax.plot(lon,lat,marker='x',markersize=20,mfc='k',mec='k')
# fudge=0.005
# ax.axis([lons.min()-fudge,lons.max()+fudge,lats.min()-fudge,lats.max()+fudge])
#
# plt.colorbar(line,ax=ax)
# plt.show(block=False)
# plt.pause(0.1)
# zs = o.history['z'][:,-2]

