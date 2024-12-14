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

data_fn = home+'data/June2024_3dhist_w_ages.p'
D = pickle.load(open(data_fn,'rb'))
for key in D.keys():
    locals()[key] = D[key]
# lon0 = -122.733651
# lat0 = 47.740037
# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773
x_center = 0.5*(x_edges[1:]+x_edges[:-1])
y_center = 0.5*(y_edges[1:]+y_edges[:-1])
z_center = 0.5*(z_edges[1:]+z_edges[:-1])

# gather mask from a random other data file, for plotting
data_fn = home+'data/f2023.06.05_test_case/ocean_his_0001.nc'
ds = nc.Dataset(data_fn)
maskr = ds['mask_rho'][:]
lonr = ds['lon_rho'][:]
latr = ds['lat_rho'][:]
xr,yr = efun.ll2xy(lonr,latr,lon0,lat0)
h= ds['h'][:]
ds.close()

# load in sampling locations
sample_info_fn = home+'data/June2024_sampling_info.csv'
df = pd.read_csv(sample_info_fn,sep=',',engine='python')

df['dt0'] = pd.to_datetime(df['date_filtered']+' '+df['time_filtered_t0']).dt.tz_localize(pdt)
df['ts0'] = df.dt0.values.astype(np.int64) // 10 ** 9
t00 = np.nanmin(df.ts0.values[:])
t11 = np.nanmax(df.ts0.values[:])
df['ts_norm'] = (df.ts0.values[:]-t00)/(t11-t00) # scale from 0 to 1
xsloc,ysloc = efun.ll2xy(df['long'].values[:],df['lat'].values[:],lon0,lat0)
df['xsloc'] = xsloc
df['ysloc'] = ysloc
xmin = np.nanmin(xsloc)
xmax = np.nanmax(xsloc)
ymin = np.nanmin(ysloc)
ymax = np.nanmax(ysloc)

Ds = {}
filter_type = df.filtered_with.unique()
for ft in filter_type:
    Ds[ft] = {}
    # I wanted to do this cleverly, but I think I need to hardcode it
    Ds[ft]['marker'] = 'o'
    Ds[ft]['markersize'] = 12
    Ds[ft]['mfc'] = 'None'
    Ds[ft]['zorder'] = 75
    if ft=='Ascension':
        Ds[ft]['marker'] = 'v'
        Ds[ft]['markersize'] = 8
        Ds[ft]['mfc'] = 'None'
        Ds[ft]['zorder']=80


ndeep = len(df[df.filtered_with=='Ascension'].sampling_location.unique()) # number of deep stations

axis_fudge=500

for t in range(len(ts_list)):
    # t=0 #while testing
    fw,fh = efun.gen_plot_props()
    fig = plt.figure(figsize=(fw*7,fh*1.25))
    # ax = fig.gca()
    gs = GridSpec(1,ndeep+1)
    axll = fig.add_subplot(gs[0])
    axs = []
    for station in range(ndeep):
        axs.append(fig.add_subplot(gs[station+1]))
    pmap = particle_map[t,-1,:,:]
    pmap[pmap==0]=np.nan
    if t==0:
        pmap_mask = ~np.isnan(pmap) # just do this once, because it shifts slightly but not in a way that matters
    axll.pcolormesh(xr,yr,maskr,cmap=cmap_mask,shading='nearest',zorder=45)
    axll.contour(xr,yr,h,levels=np.arange(0,150,10),colors='k',linewidths=0.5,zorder=50,linestyles=':')

    p=axll.pcolormesh(x_center,y_center,pmap,shading='nearest',cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(vmin=10**0,vmax=10**3))

    for ind in df.index:
        axll.plot(df.xsloc[ind],df.ysloc[ind],marker=Ds[df.filtered_with[ind]]['marker'],markersize=Ds[df.filtered_with[ind]]['markersize'],markerfacecolor='None',mec=rainbow(df.ts_norm[ind]),lw=2)
    # plt.colorbar(p)
    cbaxes = inset_axes(axll, width="4%", height="40%", loc='center right',bbox_transform=axll.transAxes,bbox_to_anchor=(-0.1,0.0,1,1))
    cbaxes.tick_params(axis='both',which='both',labelsize=8,size=2)
    cb = fig.colorbar(p, cax=cbaxes,orientation='vertical')

    axll.set_xlabel('dist from pen (m)')
    axll.set_ylabel('dist from pen (m)')
    axll.set_aspect(1)
    axll.axis([xmin-axis_fudge,xmax+axis_fudge,ymin-axis_fudge,ymax+axis_fudge])

    dt_object = datetime.fromtimestamp(ts_list[t]) # defaults to local time. Cool.
    ts_norm = (ts_list[t]-t00)/(t11-t00)
    #dt_object = datetime.utcfromtimestamp(ts_list[t]) # if you want UTC
    axll.text(0.95,0.1,dt_object.strftime("%m/%d/%Y,\n%H:%M:%S PDT"),transform=axll.transAxes,ha='right',zorder=55,fontsize=8,color=rainbow(ts_norm),fontweight='bold',bbox=dict(facecolor='white'))


    ycount =0
    yl0=0
    # for ysloc in Ds['Ascension']['ysloc']:
    for sampling_location in df[df.filtered_with=='Ascension'].sampling_location.unique():
    
        gg = df[(df['filtered_with']=='Ascension')&(df['sampling_location']==sampling_location)]

        ysloc = gg['ysloc'][:].values
        xsloc = gg['xsloc'][:].values

    
        ax= axs[ycount]
        # will want to eventually include +-50m
        # yri = np.zeros(np.shape(ysloc))
        # yci = np.zeros(np.shape(ysloc))
        # for ii in range(len(ysloc)):
        yri = np.argmin(np.abs(ysloc[-1]-yr[:,0]))
        yci = np.argmin(np.abs(ysloc[-1]-y_center))
        
        mx0 = np.argwhere(pmap_mask[yci,:])[0][0]
        mx1 = np.argwhere(pmap_mask[yci,:])[-1][0]
    
        axll.plot(x_center[mx0:mx1],[y_center[yci]]*len(x_center[mx0:mx1]),linestyle='dashed',color='k',marker='none',lw=2)
    
        pmap_z = particle_map[t,:,yci,mx0:mx1]
        ax.pcolormesh(x_center[mx0:mx1],z_center,pmap_z,shading='nearest',cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(vmin=10**0,vmax=10**3))
        ax.plot(x_center[mx0:mx1],-h[yri,mx0:mx1],linestyle='solid',color='k',marker='none')
        ax.axis([x_center[mx0],x_center[mx1],np.min(-h[yri,mx0:mx1]),3])
        yl0 = np.min([yl0,np.min(-h[yri,mx0:mx1])])
    
        ax.set_xlabel('dist from pen (m)')
        if ycount==0:
            ax.set_ylabel('depth (m)')
        ax.text(0.98,0.05,'y={:.0f} m'.format(ysloc[0]),transform=ax.transAxes,ha='right')
    
        ax.plot(x_center[mx0:mx1],[0]*len(x_center[mx0:mx1]),linestyle='dashed',color='k',marker='none',lw=2,zorder=53)
        for ind in gg.index:
            ax.plot(gg.xsloc[ind],-gg.depth_m[ind],linestyle='None',marker=Ds['Ascension']['marker'],markersize=Ds['Ascension']['markersize'],markerfacecolor='None',mec=rainbow(gg.ts_norm[ind]),zorder=500,lw=2)
    
        ycount+=1
    for ax in axs:
        ax.set_ylim([yl0,3])

    

    fig.subplots_adjust(right=0.98,left=0.05,bottom=0.11,top = 0.98,wspace=0.3)
    plt.savefig(home+f'plots/June 2024 samples/deep_time_{t}.png')
    # plt.show(block=False)
    # plt.pause(0.1)
    plt.close()


