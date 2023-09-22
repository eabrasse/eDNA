#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate a filtered time series of a lat-lon moor of ROMS output
"""

# setup
import os
import sys
pth = os.path.abspath('/Users/elizabethbrasseale/Projects/Water Quality/WQ_code/')
if pth not in sys.path:
    sys.path.append(pth)
import wqfun
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib

plt.close('all')
paired = plt.get_cmap('Paired',12)
tab10 = plt.get_cmap('tab10',10)

# load in datasets
home = '/Users/elizabethbrasseale/Projects/eDNA/'
# TG_fn = home+'data/HC_dolph_tidal_current_ATTideGauge_ks0201.p' #TG = 'tide gauge'
TG_fn = home+'data/HC_dolph_tidal_current_ATTideGauge_9445133.p'
DP_fn = home+'data/HC_dolph_tidal_current_ATBangorDolphinPen.p' #DP = 'dolphin pen'
NOAA_fn = home+'data/NOAA_MLLW_Bangor_9445133.txt' #real tide gauge data

TideGauge = pickle.load(open(TG_fn,'rb'))
DolphPen = pickle.load(open(DP_fn,'rb'))
df = pd.read_csv(NOAA_fn,engine='python',skiprows=13,delimiter=r'\t+')
df['datetime']=pd.to_datetime(df['Date ']+' '+df.Time)

# initialize figure with subplots for current and velocity data
fig = plt.figure(figsize=(14,6))
ax_ssh = fig.add_subplot(1,4,1)
ax_u = fig.add_subplot(1,4,2)
ax_v = fig.add_subplot(1,4,3)
ax_map = fig.add_subplot(1,4,4)

TideGauge['color'] = tab10(0)
NOAA_color = tab10(1)
DolphPen['color'] = tab10(2)

# mod_ls = 'solid'
# obs_ls = 'dashed'

TideGauge['ls'] = 'solid'
NOAA_ls = 'dotted'
DolphPen['ls'] = 'dashed'

TideGauge['label'] = 'Tide Gauge (model)'
NOAA_label = 'Tide Gauge (obs)'
DolphPen['label'] = 'Dolphin Pen (model)'

# first plot SSH
df['demeaned_MLLW'] = df.Pred - np.mean(df.Pred)
ax_ssh.plot(df.datetime,df.demeaned_MLLW,color=NOAA_color,linestyle=NOAA_ls,label=NOAA_label,zorder=3)

for model in TideGauge,DolphPen:
    ax_ssh.plot(model['dt_list'],model['zeta']-np.mean(model['zeta']),color=model['color'],linestyle=model['ls'],label=model['label'],zorder=1)
    ax_u.plot(model['dt_list'],model['u'],color=model['color'],linestyle=model['ls'],label=model['label'])
    ax_v.plot(model['dt_list'],model['v'],color=model['color'],linestyle=model['ls'],label=model['label'])
for ax in ax_ssh,ax_u,ax_v:
    ax.set_xlabel('Date')
    ax.set_xlim([TideGauge['dt_list'][0],TideGauge['dt_list'][-1]])
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d %H:%m"))
    plt.setp( ax.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor') 
    ax.grid()
ax_ssh.set_ylabel('SSH (m)')
ax_ssh.text(0.1,0.9,'a) Sea surface elevation (demeaned)',transform=ax_ssh.transAxes)

ax_u.set_ylabel('u (m/s)')
ax_u.text(0.1,0.9,'b) East-west velocity',transform=ax_u.transAxes)

ax_v.set_ylabel('v (m/s)')
ax_v.text(0.1,0.9,'c) North-south velocity',transform=ax_v.transAxes)

ax_ssh.legend()

# use same y lims on u and v axis
all_ylims = [ax_u.get_ylim(),ax_v.get_ylim()]
yl0 = np.min(all_ylims)
yl1 = np.max(all_ylims)
for ax in ax_u,ax_v:
    ax.set_ylim([yl0,yl1])

# add map
# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

ax_map.pcolormesh(TideGauge['lonr'],TideGauge['latr'],TideGauge['maskr'],cmap=cmap_mask)
# NOAA_lon = -(122 + 43.6/60)
# NOAA_lat = 47 + 44.9/60
NOAA_lon = -(122 + 43.8/60)
NOAA_lat = 47 + 45.5/60
markersize = 8
ax_map.plot(NOAA_lon,NOAA_lat,marker='o',mfc=NOAA_color,mec=NOAA_color,markersize=markersize+2)
lonvec = TideGauge['lonr'][0,:]
latvec = TideGauge['latr'][:,0]
for model in TideGauge,DolphPen:
    i = np.argmin(np.abs(lonvec-model['lon0']))
    j = np.argmin(np.abs(latvec-model['lat0']))
    ax_map.plot(model['lonr'][j,i],model['latr'][j,i],marker='^',mfc =model['color'],mec=model['color'],markersize=markersize)

ax_map.set_xlabel('Longitude')
ax_map.set_ylabel('Latitude')
ax_map.axis([-122.75,-122.72,47.735,47.765])
wqfun.dar(ax_map)

plt.subplots_adjust(top=0.98,left=0.05,right=0.99,wspace=0.25)
plt.show(block=False)
plt.pause(0.1)