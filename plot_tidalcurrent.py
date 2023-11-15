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
from matplotlib.gridspec import GridSpec

plt.close('all')
paired = plt.get_cmap('Paired',12)
tab10 = plt.get_cmap('tab10',10)

# load in datasets
home = '/Users/elizabethbrasseale/Projects/eDNA/'
# TG_fn = home+'data/HC_dolph_tidal_current_ATTideGauge_ks0201.p' #TG = 'tide gauge'
TG_fn = home+'data/HC_dolph_tidal_current_ATTideGauge_9445133.p'
DP_fn = home+'data/HC_dolph_tidal_current_ATBangorDolphinPen.p' #DP = 'dolphin pen'
NDP_fn = home+'data/HC_dolph_tidal_current_ATNearBangorDolphinPen.p' #DP = 'dolphin pen'
NOAA_fn = home+'data/NOAA_MLLW_Bangor_9445133.txt' #real tide gauge data

TideGauge = pickle.load(open(TG_fn,'rb'))
DolphPen = pickle.load(open(DP_fn,'rb'))
NearDolphPen = pickle.load(open(NDP_fn,'rb'))
df = pd.read_csv(NOAA_fn,engine='python',skiprows=13,delimiter=r'\t+')
df['datetime']=pd.to_datetime(df['Date ']+' '+df.Time)

# initialize figure with subplots for current and velocity data
fig = plt.figure(figsize=(14,6))
gs = GridSpec(2,4)
ax_ssh = fig.add_subplot(gs[:,0])
ax_u = fig.add_subplot(gs[0,1])
ax_v = fig.add_subplot(gs[0,2])
# ax_PCA = fig.add_subplot(gs[0,3])
ax_harm = fig.add_subplot(gs[-1,1:-1])
ax_map = fig.add_subplot(gs[:,-1])


TideGauge['color'] = tab10(0)
NOAA_color = tab10(1)
DolphPen['color'] = tab10(2)
NearDolphPen['color'] = tab10(3)

# mod_ls = 'solid'
# obs_ls = 'dashed'

TideGauge['ls'] = 'solid'
# NOAA_ls = 'dotted'
# DolphPen['ls'] = 'dashed'
# NearDolphPen['ls'] = 'dashdot'
NOAA_ls = 'solid'
DolphPen['ls'] = 'solid'
NearDolphPen['ls'] = 'solid'

TideGauge['label'] = 'Tide Gauge (model)'
NOAA_label = 'Tide Gauge (obs)'
DolphPen['label'] = 'Dolphin Pen (model)'
NearDolphPen['label'] = 'Near Dolphin Pen (model)'

lonvec = TideGauge['lonr'][0,:]
latvec = TideGauge['latr'][:,0]
TideGauge['lonnudge'] = 0
TideGauge['latnudge'] = 0
DolphPen['lonnudge'] = 0
DolphPen['latnudge'] = -2
NearDolphPen['lonnudge'] = 3
NearDolphPen['latnudge'] = -3

scale = 0.1
# first plot SSH
df['demeaned_MLLW'] = df.Pred - np.mean(df.Pred)
ax_ssh.plot(df.datetime,df.demeaned_MLLW,color=NOAA_color,linestyle=NOAA_ls,label=NOAA_label,zorder=3)
count=-1
for model in TideGauge,DolphPen,NearDolphPen:
    ax_ssh.plot(model['dt_list'],model['zeta']-np.mean(model['zeta']),color=model['color'],linestyle=model['ls'],label=model['label'],zorder=1,alpha=0.5)
    # ax_u.plot(model['dt_list'],model['u'],color=model['color'],linestyle=model['ls'],label=model['label'])
    # ax_v.plot(model['dt_list'],model['v'],color=model['color'],linestyle=model['ls'],label=model['label'])
    
    # principle component analysis
    w0 = model['u']+1j*model['v']
    w = w0 - w0.mean()
    cov = np.cov(np.imag(w),np.real(w))
    term1=cov[0,0]+cov[1,1]
    term2=np.sqrt((cov[0,0]-cov[1,1])**2 + 4*cov[1,0]**2)
    major=np.sqrt(.5*(term1+term2))
    minor=np.sqrt(.5*(term1-term2))

    theta = 0.5*np.arctan2(2*cov[1,0],(cov[0,0]-cov[1,1]))+0.5*np.pi
    model['theta']=theta
    
    w_pa = w0*np.exp(1j*theta)
    model['u_pa'] = np.real(w_pa)
    model['v_pa'] = np.imag(w_pa)
    model['u0'] = np.real(w_pa.mean())
    model['v0'] = np.imag(w_pa.mean())
    ax_u.plot(model['dt_list'],model['u_pa'],color=model['color'],linestyle=model['ls'],label=model['label'],alpha=1.,lw=0.8)
    ax_v.plot(model['dt_list'],model['v_pa'],color=model['color'],linestyle=model['ls'],label=model['label'],alpha=1.,lw=0.8)
    
    #build x axis by counting time steps in a day
    fft_frac = (model['dt_list'][1]-model['dt_list'][0])/timedelta(days=1)
    nt = len(model['dt_list'])
    nt2 = int(nt*0.5)
    freq = np.fft.fftfreq(nt,d=fft_frac)
    
    #now do fft of data
    model['fft'] = np.fft.fft(model['u_pa'])/nt
    
    #plot fft along with x-axis you built above
    ax_harm.bar(freq[:nt2]+count*0.03,model['fft'][:nt2],label=model['label'],color=model['color'],width=0.03,alpha=1)
    
    #draw ellipse
    ellx = np.linspace(-major,major,100)
    elly = np.sqrt((1-(ellx**2)/(major**2))*(minor**2))
    elly2 = -elly

    #rotate ellipse
    ellx_rot = ellx*np.cos(theta)+elly*np.sin(theta)
    elly_rot = -ellx*np.sin(theta)+elly*np.cos(theta)

    ellx2_rot = ellx*np.cos(theta)+elly2*np.sin(theta)
    elly2_rot = -ellx*np.sin(theta)+elly2*np.cos(theta)

    #plot ellipse
    i = np.argmin(np.abs(lonvec-model['lon0']))+model['lonnudge']
    j = np.argmin(np.abs(latvec-model['lat0']))+model['latnudge']
    
    ax_map.plot(model['lonr'][j,i]+scale*ellx_rot,model['latr'][j,i]+scale*elly_rot,color=model['color'],alpha=1.0,lw=0.8)
    ax_map.plot(model['lonr'][j,i]+scale*ellx2_rot,model['latr'][j,i]+scale*elly2_rot,color=model['color'],alpha=1.0,lw=0.8)
    
    count+=1
    
for ax in ax_ssh,ax_u,ax_v:
    ax.set_xlabel('Date')
    ax.set_xlim([TideGauge['dt_list'][0],TideGauge['dt_list'][-1]])
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d %H:%m"))
    plt.setp( ax.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor') 
    # ax.grid()
ax_ssh.set_ylabel('SSH (m)')
ax_ssh.text(0.1,0.9,'a) Sea surface elevation (demeaned)',transform=ax_ssh.transAxes)

ax_u.set_ylabel('u (m/s)')
# ax_u.text(0.1,0.9,'b) East-west velocity',transform=ax_u.transAxes)
ax_u.text(0.1,0.9,'b) PCA 1 velocity',transform=ax_u.transAxes)

ax_v.set_ylabel('v (m/s)')
# ax_v.text(0.1,0.9,'c) North-south velocity',transform=ax_v.transAxes)
ax_v.text(0.1,0.9,'c) PCA 2 velocity',transform=ax_v.transAxes)

ax_harm.set_xlabel('Frequency (times per day)')
ax_harm.set_xlim([-0.5,6])
ax_harm.grid()
ax_harm.text(0.1,0.9,'d) FFT breakdown of PCA1 vel (periodicity)',transform=ax_harm.transAxes)

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
markersize = 5
ax_map.plot(NOAA_lon,NOAA_lat,marker='x',mfc=NOAA_color,mec=NOAA_color,markersize=markersize+2)

for model in TideGauge,DolphPen,NearDolphPen:
    i = np.argmin(np.abs(lonvec-model['lon0']))+model['lonnudge']
    j = np.argmin(np.abs(latvec-model['lat0']))+model['latnudge']
    ax_map.plot(model['lonr'][j,i],model['latr'][j,i],marker='x',mfc =model['color'],mec=model['color'],markersize=markersize)
    # ax_map.quiver(model['lonr'][j,i],model['latr'][j,i],model['u0'],model['v0'],units='x',scale_units='x',scale=5,color=model['color'])


ax_map.set_xlabel('Longitude')
ax_map.set_ylabel('Latitude')
ax_map.axis([-122.75,-122.72,47.735,47.765])
wqfun.dar(ax_map)

plt.subplots_adjust(top=0.98,left=0.05,right=0.99,wspace=0.3,hspace=0.3)
plt.show(block=False)
plt.pause(0.1)