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
from datetime import datetime
import netCDF4 as nc
import pytz

pdt = pytz.timezone('America/Vancouver')


# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])
tab10 = plt.get_cmap('tab10',10)
Set2 = plt.get_cmap('Set2',8)
atoz = string.ascii_lowercase
rainbow = plt.get_cmap('rainbow')

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'


# load in sampling locations
drifter_fn = home+'data/microstar drifter GPS test/testmicrostarGPS.proc.txt'
drifter_df = pd.read_csv(drifter_fn,sep=',',engine='python')
phone_fn = home+'data/microstar drifter GPS test/20241010-160548 -GPSLogger Test against microstar GPS.txt'
phone_df = pd.read_csv(phone_fn,sep=',',engine='python')



fw,fh = efun.gen_plot_props()
fig = plt.figure(figsize=(fw*2,fh*1.5))

ax=fig.gca()
ax.plot(drifter_df['  Longitude'],drifter_df['   Latitude'],linestyle='none',marker='*',mec='k',mfc='none',markersize=15,label='Microstar Drifter')
ax.plot(phone_df.longitude,phone_df.latitude,linestyle='solid',lw=1,color='k',marker='None',label='Phone GPS Logger')
ax.legend()
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# fig.subplots_adjust(right=0.98,left=0.15,bottom=0.11,top = 0.98)
# plt.savefig(home+f'plots/drifter_GPS_vs_phone_GPS')
plt.show(block=False)
plt.pause(0.1)
# plt.close()


