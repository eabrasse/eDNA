import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
# from datetime import datetime, timedelta
# import pytz
from matplotlib.gridspec import GridSpec
# import matplotlib.dates as mdates
import pickle
from matplotlib import ticker
import string
import efun
import cmocean as cmo
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
from datetime import datetime, timedelta
import netCDF4 as nc
import pytz
import matplotlib.dates as mdates

tab10 = plt.get_cmap('tab10',10)
pdt = pytz.timezone('America/Vancouver')
utc = pytz.utc

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'

obs_fn = home+'data/IOOS_Wind_Bremerton9445958_20241015_20241018.csv'
dfo = pd.read_csv(obs_fn,sep=',',engine='python',skiprows=[1])
dfo['dt'] = pd.to_datetime(dfo['time'])
dfo['Wind_Direction (rad)'] = dfo.apply(lambda row: ((270-row['Wind_Direction'])*np.pi/180),axis=1)
dfo['Uwind'] = dfo.apply(lambda row: (row['Wind_Speed']*np.cos(row['Wind_Direction (rad)'])),axis=1)
dfo['Vwind'] = dfo.apply(lambda row: (row['Wind_Speed']*np.sin(row['Wind_Direction (rad)'])),axis=1)

mod_fn = home+'data/Oct2024_LOwind_Brem_Bang.p'
dmod = pickle.load(open(mod_fn,'rb'))
Brem = dmod['Brem']
Brem['col'] = tab10(2)
Bang = dmod['Bang']
Bang['col'] = tab10(1)
dmod['dt_list'] = [dt.astimezone(utc) for dt in dmod['dt_list']]

# initialize figure
fig = plt.figure(figsize=(12,6))
gs = GridSpec(2,3)
axu = fig.add_subplot(gs[0,:2])
axv = fig.add_subplot(gs[1,:2])
axmap = fig.add_subplot(gs[:,2])

# make map of stations
axmap.pcolormesh(dmod['lonr'],dmod['latr'],dmod['maskr'],cmap=cmap_mask,shading='nearest',zorder=5)
axmap.contour(dmod['lonr'],dmod['latr'],dmod['h'],levels=np.arange(0,260,20),colors='k',linewidths=0.5,zorder=6,linestyles=':')
axmap.axis([-122.874532,-122.433570,47.509171,47.858309])
efun.dar(axmap)
axmap.set_xlabel('Longitude')
axmap.set_ylabel('Latitude')
for station in Brem,Bang:
    axmap.plot(station['lon'],station['lat'],marker='o',mfc=station['col'],mec='None',markersize=10,zorder=150)
    axmap.annotate(station['name'],xy=(station['lon'],station['lat']),xytext=(station['lon']-0.02,station['lat']-0.02),\
        ha='left',va='center',fontsize=10,zorder=501)
axmap.plot(-122.626,47.562,marker='^',mfc='k',mec='None',markersize=8,zorder=160)

for ax,var in [axu,'Uwind'],[axv,'Vwind']:
    ax.plot(dfo['dt'],dfo[var],linestyle='solid',color='k')
    ax.plot(dmod['dt_list'],Brem[var],linestyle='dashed',color=Brem['col'])
    ax.plot(dmod['dt_list'],Bang[var],linestyle='dotted',color=Bang['col'])

    ax.set_xlabel('2024 Date, time (UTC)')
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d, %H:%M"))
    ax.set_ylabel('m/s')
    ax.axhline(y=0,color='gray',zorder=1)
    ax.set_xlim([dfo['dt'].min(),dfo['dt'].max()])

axu.text(0.1,0.9,'a) E-W wind speed',transform=axu.transAxes,fontsize=12)
axu.text(0.35,0.9,'Bremerton Observed',transform=axu.transAxes,color='k',fontsize=10)
axu.text(0.35,0.81,'Bremerton Modeled',transform=axu.transAxes,color=Brem['col'],fontsize=10)
axu.text(0.35,0.72,'Bangor Modeled',transform=axu.transAxes,color=Bang['col'],fontsize=10)
axv.text(0.1,0.9,'b) N-S wind speed', transform = axv.transAxes,fontsize=12)

fig.tight_layout()
plt.show(block=False)
plt.pause(0.1)