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
from scipy import interpolate

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
# home = '/Users/elizabethbrasseale/Projects/eDNA/'
home='/data2/pmr4/eab32/'


data_fn = home+'LO_data/eDNA/Oct2024_hc11_hc12_vel.p'
D = pickle.load(open(data_fn,'rb'))
for key in D.keys():
    locals()[key] = D[key]

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773


# gather mask from a random other data file, for plotting
hc = hc11
maskr = hc['maskr'][:,:]
lonr = hc['lonr'][:,:]
latr = hc['latr'][:,:]
xr,yr = efun.ll2xy(lonr,latr,lon0,lat0)
# h= hc['h'][:,:]
y0 = -1000
y1 = 500
x0 = -1500
x1 = 500
xi0 = np.argmin(np.abs(xr[0,:]-x0))
xi1 = np.argmin(np.abs(xr[0,:]-x1))
yi0 = np.argmin(np.abs(yr[:,0]-y0))
yi1 = np.argmin(np.abs(yr[:,0]-y1))
xss = 10
yss = 10

# dt_list = [utc.localize(dt) for dt in hc['dt_list']]
dt_list = hc['dt_list'][:]


surf=False
for hc in hc11,hc12:
    #stagger onto rho grid
    
    if surf: #grab -1 z index
        hc['usurf'] = 0.5*(hc['u'][:,-1,1:-1,1:]+hc['u'][:,-1,1:-1,:-1])
        hc['vsurf'] = 0.5*(hc['v'][:,-1,1:,1:-1]+hc['v'][:,-1,:-1,1:-1])
        plottitle = 'Surface velocity'
    else: # use u/v bar
        hc['usurf'] = 0.5*(hc['ubar'][:,1:-1,1:]+hc['ubar'][:,1:-1,:-1])
        hc['vsurf'] = 0.5*(hc['vbar'][:,1:,1:-1]+hc['vbar'][:,:-1,1:-1])
        plottitle = 'Vertically-averaged velocity'
    
    hc['usurf'][hc['usurf']>50] = np.nan
    hc['vsurf'][hc['vsurf']>50] = np.nan

    hc['speed'] = np.sqrt(hc['usurf']**2+hc['vsurf']**2)
    hc['speed'][hc['speed']>500]= 0

    hc['cmap'] = cmo.cm.speed
    hc['vmin'] = 0
    hc['vmax'] = 1
    
hcdiff = hc12.copy() # use masked delta pier grid for the difference
hcdiff['grid'] = 'difference'
vname_list = ['speed','usurf','vsurf']
for vname in vname_list:
    hcdiff[vname] = hc12[vname]-hc11[vname]
hcdiff['cmap'] = cmo.cm.balance
hcdiff['vmin'] = -0.25
hcdiff['vmax'] = 0.25

hc_list = [hc11,hc12,hcdiff]


    
for t in range(len(hc['dt_list'])):
    # t=15 #while testing
    fw,fh = efun.gen_plot_props()
    fig = plt.figure(figsize=(fw*1.4*len(hc_list),fh*1.1))
    gs = GridSpec(1,len(hc_list))
    count=0
    for hc in hc_list:
        axll = fig.add_subplot(gs[count])
        # axll = fig.gca()
        axll.set_aspect(1)
        axll.axis([-1500,500,-1000,500])

        axll.pcolormesh(xr,yr,hc['maskr'],cmap=cmap_mask,shading='nearest',zorder=5)


        p=axll.pcolormesh(xr[1:-1,1:-1],yr[1:-1,1:-1],hc['speed'][t,:,:],shading='nearest',cmap=hc['cmap'],vmin=hc['vmin'],vmax=hc['vmax'],alpha=0.75)
        axll.quiver(xr[yi0:yi1:yss,xi0:xi1:xss],yr[yi0:yi1:yss,xi0:xi1:xss],hc['usurf'][t,yi0:yi1:yss,xi0:xi1:xss],hc['vsurf'][t,yi0:yi1:yss,xi0:xi1:xss],
            color='k',scale_units='xy',scale=0.002)


        if count!=1:
            cbaxes = inset_axes(axll, width="4%", height="40%", loc='center right',bbox_transform=axll.transAxes,bbox_to_anchor=(-0.12,-0.08,1,1))
            cbaxes.tick_params(axis='both',which='both',labelsize=8,size=2)
            cb = fig.colorbar(p, cax=cbaxes,orientation='vertical')
            cbaxes.set_ylabel('Speed (m/s)')
        
            axll.quiver(-500,-800,0.1,0.1,color='k',scale_units='xy',scale=0.002,zorder=50)
            axll.text(-500,-900,'0.1 m/s',zorder=600)

        axll.set_xlabel('Dist from pen (m)')



        # dt_object = datetime.fromtimestamp(ts_list[t]) # defaults to local time. Cool.
        axll.text(0.1,0.9,hc['grid'],transform=axll.transAxes,ha='left',fontsize=12,zorder=600,bbox=dict(facecolor='white'))

        if count==0:
            axll.text(0.1,1.05,plottitle,transform=axll.transAxes,ha='left',fontsize=12,zorder=600)
            axll.text(0.95,0.1,dt_list[t].astimezone(pdt).strftime("%m/%d/%Y,\n%H:%M PDT"),transform=axll.transAxes,ha='right',zorder=55,fontsize=8,color='k',fontweight='bold',bbox=dict(facecolor='white'))
            axll.set_ylabel('Dist from pen (m)')
        else:
            axll.set_yticklabels([''])
        count+=1


    fig.subplots_adjust(right=0.98,left=0.1,bottom=0.10,top = 0.95,wspace=0.1)
    plt.savefig(home+f'etools/plots/delta pier mask/problem_vel_bar_{t:0>2}.png')
    # plt.show(block=False)
    # plt.pause(0.1)
    plt.close()


