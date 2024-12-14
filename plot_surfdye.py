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

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

atoz = string.ascii_lowercase

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'

data_fn = home+'data/HC_surfdye_0605.p'

lon0 = -122.733651
lat0 = 47.740037

fw,fh = efun.gen_plot_props()
fig = plt.figure(figsize=(fw*2,fh))
gs = GridSpec(1,2)

Dp = pickle.load(open(data_fn,'rb'))

# gather some fields, for convenience
lonr = Dp['lonr'][:]
latr = Dp['latr'][:]
maskr = Dp['maskr'][:]

dye = np.ma.masked_where(Dp['dye']>10,Dp['dye'])

Dmax = {}
Dmax['var'] = np.max(dye,axis=0)
Dmax['ax'] = plt.subplot(gs[0])
Dmax['norm'] = matplotlib.colors.LogNorm(vmax=1e-1,vmin=1e-5)
Dmax['levels'] = [1e-3]
Dmax['label'] = 'Max dye concentration'
Dmax['cmap'] = cmo.cm.matter

Dtime = {}
perctime= np.sum(dye>1e-3,axis=0)/24
perctime[perctime==0] = np.nan
Dtime['var'] =perctime
Dtime['ax'] = plt.subplot(gs[1])
Dtime['norm'] = matplotlib.colors.Normalize(vmin=0, vmax=1)
Dtime['label'] = r'$\%$ time > 1e-3'
Dtime['levels'] = [0.25,0.5,0.75]
Dtime['cmap'] = cmo.cm.dense

xgrid,ygrid = efun.ll2xy(lonr,latr,lon0,lat0)
nudge = 150

count=0
for D in Dmax,Dtime:
    D['ax'].pcolormesh(xgrid,ygrid,maskr,cmap=cmap_mask,shading='nearest',zorder=100)
    p=D['ax'].pcolormesh(xgrid,ygrid,D['var'],cmap=D['cmap'],norm=D['norm'],zorder=1,shading='nearest')
    
    D['ax'].contour(xgrid,ygrid,maskr,levels=[0.5],colors=['k'],linewidths=[1.5],zorder=150)
    D['ax'].set_xlabel('Dist from floating pen (m)')
    D['ax'].set_ylabel('Dist from floating pen (m)')
    cbaxes = inset_axes(D['ax'], width="4%", height="40%", loc='lower right',bbox_transform=D['ax'].transAxes,bbox_to_anchor=(-0.28,0,1,1))
    D['cb'] = fig.colorbar(p, cax=cbaxes, orientation='vertical')
    
    pl=D['ax'].contour(xgrid,ygrid,D['var'],colors=['k'],linewidths=[0.5],levels=D['levels'],zorder=2)
    D['cb'].add_lines(pl)
    
    D['ax'].set_aspect(1)
    D['ax'].axis([-nudge,nudge,-nudge,nudge])
    D['ax'].plot(0,0,marker='x',color='yellow')
    D['ax'].text(0.1,0.9,atoz[count]+') '+ D['label'],transform=D['ax'].transAxes,ha='left',va='top',zorder=5000)
    
    # D['ax'].grid(zorder=-1)
    
    count+=1


# pl=Dmax['ax'].contour(xgrid,ygrid,Dmax['var'],colors=['k'],linewidths=[0.5],levels=Dmax['levels'],zorder=2)
# Dmax['cb'].add_lines(pl)

fig.subplots_adjust(right=0.98,left=0.11,bottom=0.15,top = 0.98,wspace=0.3)

plt.show(block=False)
plt.pause(0.1)
# plt.close()


