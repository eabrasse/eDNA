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
utc = pytz.utc


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


data_fn = home+'/data/Oct2024_1845rel_surfdye_15min.p'
D = pickle.load(open(data_fn,'rb'))
for key in D.keys():
    locals()[key] = D[key]
# lon0 = -122.733651
# lat0 = 47.740037
# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773
# surface_dye = np.ma.masked_where(surface_dye,surface_dye>1000)

# dt = ts_list[1]-ts_list[0]
# particle_map = particle_map/dv

# gather mask from a random other data file, for plotting
data_fn = home+'data/Oct2024_HC2wind.p'
ds = nc.Dataset(data_fn)
maskr = ds['mask_rho'][1:-1,1:-1]
lonr = ds['lon_rho'][1:-1,1:-1]
latr = ds['lat_rho'][1:-1,1:-1]
xr,yr = efun.ll2xy(lonr,latr,lon0,lat0)
h= ds['h'][1:-1,1:-1]
ds.close()


for t in range(len(dt_list)):
    # t=15 #while testing
    fw,fh = efun.gen_plot_props()
    fig = plt.figure(figsize=(fw*(2+nstation),fh*1))
    # gs = GridSpec(1,2+nstation)
    # axll = fig.add_subplot(gs[:2])
    axll = fig.gca()


    axll.pcolormesh(xr,yr,maskr,cmap=cmap_mask,shading='nearest',zorder=5)
    axll.contour(xr,yr,h,levels=np.arange(0,150,10),colors='k',linewidths=0.5,zorder=6,linestyles=':')

    p=axll.pcolormesh(xr,yr,surface_dye[t,1:-1,1:-1],shading='nearest',cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(vmin=1e-5,vmax=0.05))

    cbaxes = inset_axes(axll, width="4%", height="40%", loc='center right',bbox_transform=axll.transAxes,bbox_to_anchor=(-0.1,-0.05,1,1))
    cbaxes.tick_params(axis='both',which='both',labelsize=8,size=2)
    cb = fig.colorbar(p, cax=cbaxes,orientation='vertical')

    axll.set_xlabel('Dist from pen (m)')
    axll.set_ylabel('Dist from pen (m)')
    axll.set_aspect(1)
    axll.axis([-1500,500,-1000,500])

    # dt_object = datetime.fromtimestamp(ts_list[t]) # defaults to local time. Cool.
    axll.text(0.95,0.1,dt_list[t].astimezone(pdt).strftime("%m/%d/%Y,\n%H:%M PDT"),transform=axll.transAxes,ha='right',zorder=55,fontsize=8,color='k',fontweight='bold',bbox=dict(facecolor='white'))
    axll.text(0.1,0.9,'Dye',transform=axll.transAxes,ha='left',fontsize=12,zorder=600)


    fig.subplots_adjust(right=0.98,left=0.1,bottom=0.15,top = 0.95,wspace=0.5)
    plt.savefig(home+f'plots/Oct 2024 samples/dye 1845 rel/dye_station_profiles_{t:0>2}.png')
    # plt.show(block=False)
    # plt.pause(0.1)
    plt.close()


