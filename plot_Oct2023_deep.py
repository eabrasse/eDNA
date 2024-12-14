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

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'

data_fn = home+'data/Oct2023_3dhist_w_ages_VL_only.p'
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
Ds = {}
filter_type = df.filtered_with.unique()
# df['dt0'] = pd.to_datetime(df[['date_filtered','time_filtered_0']])
# df['ts0'] = df.dt0.values.astype(np.int64) // 10 ** 9

xmin=0
xmax=0
ymin=0
ymax=0

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
    
    gg = df[df['filtered_with']==ft]
    Ds[ft]['sampling_locations'] = gg.sampling_location.unique()

    for sampling_location in Ds[ft]['sampling_locations']:
        ggg = gg[gg['sampling_location']==sampling_location]
        Ds[ft][sampling_location] = {}
        # all data entries at each sampling location should have the same lat/lon
        # so just read the first index
        lon = ggg['long'].values[:]
        lat = ggg['lat'].values[:]
        
        Ds[ft][sampling_location]['xsloc'],Ds[ft][sampling_location]['ysloc'] = efun.ll2xy(lon,lat,lon0,lat0)
        xmin = np.nanmin([xmin,np.nanmin(Ds[ft][sampling_location]['xsloc'])])
        xmax = np.nanmax([xmax,np.nanmax(Ds[ft][sampling_location]['xsloc'])])
        ymin = np.nanmin([ymin,np.nanmin(Ds[ft][sampling_location]['ysloc'])])
        ymax = np.nanmax([ymax,np.nanmax(Ds[ft][sampling_location]['ysloc'])])

        # if ft=='Ascension':
        #     hh = df[(df['filtered_with']==ft)&(df['sampling_location']==sampling_location)]
        #     depths = hh.depth_m.unique()
        #     xdepths = np.zeros(depths.shape)
        #     ydepths = np.zeros(depths.shape)
        #     dcount=0
        #     for depth in depths:
        #         lon_depth = hh[hh['depth_m']==depth]['long'].values[0]
        #         lat_depth = hh[hh['depth_m']==depth]['lat'].values[0]
        #         xdepths[dcount],ydepths[dcount] = efun.ll2xy(lon_depth,lat_depth,lon0,lat0)
        #         dcount+=1


ndeep = len(Ds['Ascension']['sampling_locations']) # number of deep stations

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

    # for ft in filter_type:
        # for sampling_location in Ds[ft]['sampling_locations']:
            # axll.plot(Ds[ft][sampling_location]['xsloc'][:],Ds[ft][sampling_location]['ysloc'][:],linestyle='None',marker=Ds[ft]['marker'],markersize=Ds[ft]['markersize'],markerfacecolor=Ds[ft]['mfc'],mec='k',zorder=Ds[ft]['zorder'])
    # plt.colorbar(p)
    cbaxes = inset_axes(axll, width="4%", height="40%", loc='center right',bbox_transform=axll.transAxes,bbox_to_anchor=(-0.1,0.0,1,1))
    cbaxes.tick_params(axis='both',which='both',labelsize=8,size=2)
    cb = fig.colorbar(p, cax=cbaxes,orientation='vertical')

    axll.set_xlabel('dist from pen (m)')
    axll.set_ylabel('dist from pen (m)')
    axll.set_aspect(1)
    axll.axis([-1500,500,-1500,1000])

    dt_object = datetime.fromtimestamp(ts_list[t]) # defaults to local time. Cool.
    #dt_object = datetime.utcfromtimestamp(ts_list[t]) # if you want UTC
    axll.text(0.95,0.1,dt_object.strftime("%m/%d/%Y,\n%H:%M:%S PDT"),transform=axll.transAxes,ha='right',zorder=55,fontsize=8)


    ycount =0
    yl0=0
    # for ysloc in Ds['Ascension']['ysloc']:
    for ysloc in [-1000,-500,0]:
    
        # gg = df[(df['filtered_with']=='Ascension')&(df['sampling_location']==sampling_location)]
        # depths = gg.depth_m.unique()
        # xdepths = np.zeros(depths.shape)
        # ydepths = np.zeros(depths.shape)
        # dcount=0
        # for depth in depths:
        #     lon_depth = gg[gg['depth_m']==depth]['long'].values[0]
        #     lat_depth = gg[gg['depth_m']==depth]['lat'].values[0]
        #     xdepths[dcount],ydepths[dcount] = efun.ll2xy(lon_depth,lat_depth,lon0,lat0)
        #     dcount+=1
        # ydiff = ydepths.max()-ydepths.min()
        # print(f'y spread ={ydiff:} m')
    
        # ysloc = Ds['Ascension'][sampling_location]['ysloc'][:]
        # xsloc = Ds['Ascension'][sampling_location]['xsloc'][:]

    
        ax= axs[ycount]
        # will want to eventually include +-50m
        # yri = np.zeros(np.shape(ysloc))
        # yci = np.zeros(np.shape(ysloc))
        # for ii in range(len(ysloc)):
        yri = np.argmin(np.abs(ysloc-yr[:,0]))
        yci = np.argmin(np.abs(ysloc-y_center))
        
        mx0 = np.argwhere(maskr[yri,:])[0][0]
        mx1 = np.argwhere(maskr[yri,:])[-1][0]
    
        axll.plot(x_center[mx0:mx1],[y_center[yci]]*len(x_center[mx0:mx1]),linestyle='dashed',color=Set2(ycount),marker='none',lw=2)
    
        pmap_z = particle_map[t,:,yci,mx0:mx1]
        ax.pcolormesh(x_center[mx0:mx1],z_center,pmap_z,shading='nearest',cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(vmin=10**0,vmax=10**3))
        ax.plot(x_center[mx0:mx1],-h[yri,mx0:mx1],linestyle='solid',color='k',marker='none')
        ax.axis([x_center[mx0],x_center[mx1],np.min(-h[yri,mx0:mx1]),3])
        yl0 = np.min([yl0,np.min(-h[yri,mx0:mx1])])
    
        ax.set_xlabel('dist from pen (m)')
        if ycount==0:
            ax.set_ylabel('depth (m)')
        ax.text(0.98,0.05,'y={:.0f} m'.format(ysloc),transform=ax.transAxes,ha='right')
    
        ax.plot(x_center[mx0:mx1],[0]*len(x_center[mx0:mx1]),linestyle='dashed',color=Set2(ycount),marker='none',lw=2,zorder=53)
        # ax.plot(xdepths,-depths,linestyle='None',marker=Ds['Ascension']['marker'],markersize=Ds['Ascension']['markersize'],markerfacecolor=Ds['Ascension']['mfc'],mec='k',zorder=500)
    
        ycount+=1
    for ax in axs:
        ax.set_ylim([yl0,3])

    

    fig.subplots_adjust(right=0.98,left=0.05,bottom=0.11,top = 0.98,wspace=0.3)
    plt.savefig(home+f'plots/Oct 2023 planning/deep_VL_only_{t}.png')
    # plt.show(block=False)
    # plt.pause(0.1)
    plt.close()


