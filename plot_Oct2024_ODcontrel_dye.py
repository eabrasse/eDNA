import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
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

ncfile = home+'data/opendrift/HCdolph_cont_rel_test_hdiff02.nc'
dsod = nc.Dataset(ncfile)
ot = dsod['time'][:]
ot_list = []
for ott in ot:
    ot_list.append(datetime(1970,1,1,tzinfo=utc)+timedelta(seconds=ott))
 

data_fn = home+'data/Oct2024_1845rel_surfdye_15min.p'
D = pickle.load(open(data_fn,'rb'))
for key in D.keys():
    locals()[key] = D[key]
# lon0 = -122.733651
# lat0 = 47.740037
# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773
# surface_dye = np.ma.masked_where(surface_dye,surface_dye>1000)
dt_list = D['dt_list'][:]

shared_dt = [dt for dt in dt_list if dt in ot_list]


# gather mask from a random other data file, for plotting
data_fn = home+'data/f2023.06.05_test_case/ocean_his_0001.nc'
ds = nc.Dataset(data_fn)
maskr = ds['mask_rho'][1:-1,1:-1]
lonr = ds['lon_rho'][1:-1,1:-1]
latr = ds['lat_rho'][1:-1,1:-1]
lonp = ds['lon_psi'][1:-1,1:-1]
latp = ds['lat_psi'][1:-1,1:-1]
xr,yr = efun.ll2xy(lonr,latr,lon0,lat0)
x_edges,y_edges = efun.ll2xy(lonp[0,:],latp[:,0],lon0,lat0)
h= ds['h'][1:-1,1:-1]
ds.close()

# for plotting: import fence coordinates
fence_fn = home+'data/perimeter_fence_coords.csv'
fence_df = pd.read_csv(fence_fn,sep=',',engine='python')
fence_x,fence_y = efun.ll2xy(fence_df.lon_decimal,fence_df.lat_decimal,lon0,lat0)

xfence_fn = home+'data/perimeter_fence_coords_X.csv'
xfence_df = pd.read_csv(xfence_fn,sep=',',engine='python')
xfence_x,xfence_y = efun.ll2xy(xfence_df.lon_decimal,xfence_df.lat_decimal,lon0,lat0)

px,py = efun.ll2xy(dsod['lon'][:],dsod['lat'][:],lon0,lat0)
nstation = len(station_list)

# dv_station = dz*np.pi*(100**2)
# for station in station_list:
    # station['profile'][station['profile']==0] = 1e-10

[VL,HA,NB] = station_list
# # VL pen and Husbandry area sampling stations
# VL = {}
# VL['name'] = 'VL\nPen\nSouth'
# VL['lon'] = -122.733598
# VL['lat'] = 47.740000
# VL['x'], VL['y'] = efun.ll2xy(VL['lon'],VL['lat'],lon0,lat0)
# # VL_particle_profile = np.zeros((nt,nz-2))
# VL['xi'] = np.argmin(np.abs(xr[0,1:-1]-VL['x']))
# VL['yi'] = np.argmin(np.abs(yr[1:-1,0]-VL['y']))
VL['col'] = tab10(8)
VL['text_x'] = 125
VL['text_y'] = 0
# VL['profile'] = VL_particle_profile
#
# HA = {}
# HA['name'] = 'Husbandry\nArea'
# HA['lat'] = 47.742236
# HA['lon'] = -122.729975
# HA['x'], HA['y'] = efun.ll2xy(HA['lon'],HA['lat'],lon0,lat0)
# # HA_particle_profile = np.zeros((nt,nz-2))
# HA['xi'] = np.argmin(np.abs(xr[0,1:-1]-HA['x']))
# HA['yi'] = np.argmin(np.abs(yr[1:-1,0]-HA['y']))
HA['col'] = tab10(6)
HA['text_x'] = 0
HA['text_y'] = 125
# HA['profile'] = HA_particle_profile
#
# NB = {}
# NB['name'] = 'NOAA\nBoat'
# NB['lat'] = 47.736613
# NB['lon'] = -122.743109
# NB['x'],NB['y'] = efun.ll2xy(NB['lon'],NB['lat'],lon0,lat0)
# NB['xi'] = np.argmin(np.abs(xr[0,1:-1]-NB['x']))
NB['col'] = tab10(9)
NB['text_x'] = -150
NB['text_y'] = 0
# NB['profile'] = NB_particle_profile

# for t in range(len(dt_list)):
count=0
for dtt in shared_dt:
    # dtt = shared_dt[15]
    t = np.argmin([np.abs(dt-dtt) for dt in dt_list])
    to = np.argmin([np.abs(ot-dtt) for ot in ot_list])

    fw,fh = efun.gen_plot_props()
    fig = plt.figure(figsize=(fw*3,fh*1))
    # gs = GridSpec(1,2+nstation)
    gs = GridSpec(1,2)
    axll = fig.add_subplot(gs[0])
    axod = fig.add_subplot(gs[1])
    # axll = fig.gca()
    # for sta in range(nstation):
        # station_list[sta]['ax'] = fig.add_subplot(gs[2+sta])

    # zmask = z_center[t,:,:,:]>-10
    # pmap = np.sum(zmask*particle_map[t,:,:,:],axis=0)/np.sum(zmask*dv[t,:,:,:],axis=0)
    # pmap[pmap==0]=np.nan
    for ax in axll, axod:
        ax.pcolormesh(xr,yr,maskr,cmap=cmap_mask,shading='nearest',zorder=5)
        ax.contour(xr,yr,h,levels=np.arange(0,150,10),colors='k',linewidths=0.5,zorder=6,linestyles=':')
        ax.set_xlabel('Dist from pen (m)')
        ax.set_ylabel('Dist from pen (m)')
        ax.set_aspect(1)
        ax.axis([-1500,500,-1000,500])
        ax.plot(fence_x,fence_y,linestyle='dashed',marker='o',color='k',mfc='k',mec='k',markersize=2,lw=1)
        ax.plot(xfence_x,xfence_y,linestyle='dashed',marker='o',color='k',mfc='k',mec='k',markersize=2,lw=1)
        ax.plot(VL['x'],VL['y'],marker='o',linestyle='none',mec='k',mfc=VL['col'],markersize=8,zorder=70)
    axll.text(VL['x']+VL['text_x'],VL['y']+VL['text_y'],VL['name'],fontweight='bold',color='k',zorder=550,ha='center',va='center')

    p=axll.pcolormesh(xr,yr,surface_dye[t,1:-1,1:-1],shading='nearest',cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(vmin=1e-5,vmax=0.05))

    cbaxes = inset_axes(axll, width="4%", height="40%", loc='center right',bbox_transform=axll.transAxes,bbox_to_anchor=(-0.1,-0.05,1,1))
    cbaxes.tick_params(axis='both',which='both',labelsize=8,size=2)
    cb = fig.colorbar(p, cax=cbaxes,orientation='vertical')

    hist, xedges, yedges = np.histogram2d(px[:,to],py[:,to],bins=[x_edges,y_edges])
    p=axod.pcolormesh(xr[1:-1,1:-1],yr[1:-1,1:-1],hist.T,shading='nearest',cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(vmin=1,vmax=100))

    cbaxes = inset_axes(axod, width="4%", height="40%", loc='center right',bbox_transform=axod.transAxes,bbox_to_anchor=(-0.1,-0.05,1,1))
    cbaxes.tick_params(axis='both',which='both',labelsize=8,size=2)
    cb = fig.colorbar(p, cax=cbaxes,orientation='vertical')

    # add particles
    # axll.scatter(px[:,to],py[:,to],c='k',s=5,alpha=1,zorder=1000)


    # dt_object = datetime.fromtimestamp(ts_list[t]) # defaults to local time. Cool.
    axll.text(0.95,0.1,dt_list[t].astimezone(pdt).strftime("%m/%d/%Y,\n%H:%M PDT"),transform=axll.transAxes,ha='right',zorder=55,fontsize=8,color='k',fontweight='bold',bbox=dict(facecolor='white'))
    axll.text(0.1,0.9,'Dye',transform=axll.transAxes,ha='left',fontsize=12,zorder=600)
    axod.text(0.1,0.9,'Opendrift',transform=axod.transAxes,ha='left',fontsize=12,zorder=600)




    #add icons at locations
    # for station in station_list:
    #     axll.plot(station['x'],station['y'],marker='o',linestyle='none',mec='k',mfc=station['col'],markersize=8,zorder=70)
    #     axll.text(station['x']+station['text_x'],station['y']+station['text_y'],station['name'],fontweight='bold',color='k',zorder=550,ha='center',va='center')
    #     station['ax'].plot(station['profile'][t,:]+1e-10,station['zr'][t,:],linestyle='solid',marker='none',color=station['col'],lw=2)
    #     station['ax'].set_xlabel('Dye')
    #     station['ax'].set_ylabel('Depth (m)')
    #     station['ax'].set_ylim(-50,2)
    #     station['ax'].set_xlim(1e-11,1e-2)
    #     station['ax'].set_xscale('log')
    #     # xticks = station['ax'].get_xticks()
    #     # xticklabels = station['ax'].get_xticklabels()
    #     # xticklabels[0] = '0'
    #     station['ax'].set_xticks([1e-10,1e-8,1e-6,1e-4,1e-2],['0','$\\mathdefault{10^{-8}}$','$\\mathdefault{10^{-6}}$','$\\mathdefault{10^{-4}}$','$\\mathdefault{10^{-2}}$'])
    #
    #     station['ax'].spines[['right', 'top']].set_visible(False)
    #     station['ax'].text(0.9,0.8,station['name'],transform=station['ax'].transAxes,ha='right',va='center')

    fig.subplots_adjust(right=0.98,left=0.1,bottom=0.15,top = 0.95,wspace=0.5)
    plt.savefig(home+f'plots/Oct 2024 samples/dye vs opendrift/fig_hdiff02_2panel{count:0>2}.png')
    # plt.show(block=False)
    # plt.pause(0.1)
    plt.close()
    count+=1


