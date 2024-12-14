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

pdt = pytz.timezone('America/Vancouver')
utc = pytz.utc


# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])
tab10 = plt.get_cmap('tab10',10)
Set2 = plt.get_cmap('Set2',8)
tab20b = plt.get_cmap('tab20b',20)
atoz = string.ascii_lowercase
rainbow = plt.get_cmap('rainbow')

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

dye_0lag = {}
dye_1lag = {}
particles = {}


dye_0lag['data_fn'] = home+'/data/Oct2024_1845rel_surfdye_15min.p'
dye_0lag['label'] = 'Continuous dye release'
dye_0lag['ls'] = 'solid'
dye_1lag['data_fn'] = home+'/data/Oct2024_1945rel_surfdye_15min.p'
dye_1lag['label'] = 'Continuous dye release\n1 hr lag'
dye_1lag['ls'] = 'dashed'
particles['data_fn'] = home+'/data/Oct2024_3dhist_zw_no_ages_15min_VL_only.p'
particles['label'] = 'Hourly particle release'
particles['ls'] = 'dotted'

model_list = [dye_0lag,dye_1lag,particles]

for mod in model_list:
    D = pickle.load(open(mod['data_fn'],'rb'))
    for key in D.keys():
        mod[key] = D[key]

# create comparable time stamps
particles['dt_list'] = [utc.localize(datetime.utcfromtimestamp(ts)) for ts in particles['ts_list']]
dye_1lag['dt_list0'] = dye_1lag['dt_list'][:]
dye_1lag['dt_list'] = [dt-timedelta(hours=1) for dt in dye_1lag['dt_list0']]

#make master time list containing only elements common to all 3 models

dt0_list = [min(mod['dt_list']) for mod in model_list]
dt0 = max(dt0_list)
dt1_list = [max(mod['dt_list']) for mod in model_list]
dt1 = min(dt1_list)

dt_list = pd.date_range(start=dt0,end = dt1, freq="15min").to_pydatetime().tolist()

#some grid info for particles plot
x_center = 0.5*(particles['x_edges'][1:]+particles['x_edges'][:-1])
y_center = 0.5*(particles['y_edges'][1:]+particles['y_edges'][:-1])
z_center = 0.5*(particles['z_edges_15min'][:,1:,:,:]+particles['z_edges_15min'][:,:-1,:,:])
z_center = 0.5*(z_center[:,:,:,1:]+z_center[:,:,:,:-1])
z_center = 0.5*(z_center[:,:,1:,:]+z_center[:,:,:-1,:])
xx,yy = np.meshgrid(x_center,y_center)
dx = np.diff(particles['x_edges']) #not uniform
dy = np.diff(particles['y_edges']) # not uniform

dz = np.diff(particles['z_edges_15min'],axis=1)
dz = 0.5*(dz[:,:,:,1:]+dz[:,:,:,:-1])
dz = 0.5*(dz[:,:,1:,:]+dz[:,:,:-1,:])
nt,nz,ny,nx = np.shape(dz)
dx = np.tile(np.reshape(dx,[1,1,1,nx]),[nt,nz,ny,1])
dy = np.tile(np.reshape(dy,[1,1,ny,1]),[nt,nz,1,nx])
dv = dx*dy*dz #4d array


dv_station = dz*np.pi*(100**2)
for station in particles['station_list']:
    station['profile'] = station['profile'][:,:]/dv_station[:,:,station['yi'],station['xi']]
    station['zr'] = z_center[:,:,station['yi'],station['xi']]


# gather mask from a random other data file, for plotting
data_fn = home+'data/f2023.06.05_test_case/ocean_his_0001.nc'
ds = nc.Dataset(data_fn)
maskr = ds['mask_rho'][1:-1,1:-1]
lonr = ds['lon_rho'][1:-1,1:-1]
latr = ds['lat_rho'][1:-1,1:-1]
xr,yr = efun.ll2xy(lonr,latr,lon0,lat0)
h= ds['h'][1:-1,1:-1]
ds.close()

# for plotting: import fence coordinates
fence_fn = home+'data/perimeter_fence_coords.csv'
fence_df = pd.read_csv(fence_fn,sep=',',engine='python')
fence_x,fence_y = efun.ll2xy(fence_df.lon_decimal,fence_df.lat_decimal,lon0,lat0)

xfence_fn = home+'data/perimeter_fence_coords_X.csv'
xfence_df = pd.read_csv(xfence_fn,sep=',',engine='python')
xfence_x,xfence_y = efun.ll2xy(xfence_df.lon_decimal,xfence_df.lat_decimal,lon0,lat0)

nstation = len(particles['station_list'])


[VL,HA,NB] = particles['station_list']

VL['col'] = 0
VL['text_x'] = 125
VL['text_y'] = 0

HA['col'] = 4
HA['text_x'] = 0
HA['text_y'] = 125

NB['col'] = 12
NB['text_x'] = -150
NB['text_y'] = 0


for t in range(len(dt_list)):
    # t=15 #while testing
    fw,fh = efun.gen_plot_props()
    fig = plt.figure(figsize=(fw*4,fh*2))
    dye_0lag['ax'] = fig.add_subplot(2,3,1)
    dye_1lag['ax'] = fig.add_subplot(2,3,2)
    particles['ax'] = fig.add_subplot(2,3,3)
    VL['ax'] = fig.add_subplot(2,3,4)
    HA['ax'] = fig.add_subplot(2,3,5)
    NB['ax'] = fig.add_subplot(2,3,6)

    count=0

    for mod in model_list:
        tind = mod['dt_list'].index(dt_list[t])

        mod['ax'].pcolormesh(xr,yr,maskr,cmap=cmap_mask,shading='nearest',zorder=5)
        mod['ax'].contour(xr,yr,h,levels=np.arange(0,150,10),colors='k',linewidths=0.5,zorder=6,linestyles=':')

        if 'surface_dye' in mod.keys(): # dye
            p=mod['ax'].pcolormesh(xr,yr,mod['surface_dye'][tind,1:-1,1:-1],shading='nearest',cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(vmin=1e-5,vmax=0.05))
    
        if 'particle_map' in mod.keys(): #particles

            zmask = z_center[t,:,:,:]>-10
            pmap = np.sum(zmask*mod['particle_map'][tind,:,:,:],axis=0)/np.sum(zmask*dv[tind,:,:,:],axis=0)
            pmap[pmap==0]=np.nan

            p=mod['ax'].pcolormesh(x_center,y_center,pmap,shading='nearest',cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(vmin=10**-3,vmax=3*10**0))

        cbaxes = inset_axes(mod['ax'], width="4%", height="40%", loc='center right',bbox_transform=mod['ax'].transAxes,bbox_to_anchor=(-0.1,-0.05,1,1))
        cbaxes.tick_params(axis='both',which='both',labelsize=8,size=2)
        cb = fig.colorbar(p, cax=cbaxes,orientation='vertical')

        mod['ax'].set_xlabel('Dist from pen (m)')
        mod['ax'].set_ylabel('Dist from pen (m)')
        mod['ax'].set_aspect(1)
        mod['ax'].axis([-1500,500,-1000,500])

        mod['ax'].plot(fence_x,fence_y,linestyle='dashed',marker='o',color='k',mfc='k',mec='k',markersize=2,lw=1)
        mod['ax'].plot(xfence_x,xfence_y,linestyle='dashed',marker='o',color='k',mfc='k',mec='k',markersize=2,lw=1)

        for stat in range(nstation):
            station = mod['station_list'][stat]
            pstation = particles['station_list'][stat]
        
            mod['ax'].plot(station['x'],station['y'],marker='o',linestyle='none',mec=tab20b(pstation['col']),mfc=tab20b(pstation['col']+1),markersize=8,zorder=70)
        
            pstation['ax'].plot(station['profile'][t,:]+1e-10,station['zr'][t,:],linestyle=mod['ls'],marker='none',color=tab20b(pstation['col']+count),lw=2,label=mod['label'])
            pstation['ax'].set_xlabel('Concentration')
            pstation['ax'].set_ylabel('Depth (m)')
            pstation['ax'].set_ylim(-50,2)
            pstation['ax'].set_xlim(1e-11,1e-2)
            pstation['ax'].set_xscale('log')
            pstation['ax'].set_xticks([1e-10,1e-8,1e-6,1e-4,1e-2],['0','$\\mathdefault{10^{-8}}$','$\\mathdefault{10^{-6}}$','$\\mathdefault{10^{-4}}$','$\\mathdefault{10^{-2}}$'])

            pstation['ax'].spines[['right', 'top']].set_visible(False)
            pstation['ax'].text(0.9,0.8,station['name'],transform=pstation['ax'].transAxes,ha='right',va='center')

        mod['ax'].text(0.1,0.9,mod['label'],transform=mod['ax'].transAxes,ha='left',va='top',fontsize=12,zorder=600)
    
        count+=1

    # add annotation only to first axis 
    for station in particles['station_list']:      
        dye_0lag['ax'].text(station['x']+station['text_x'],station['y']+station['text_y'],station['name'],fontweight='bold',color='k',zorder=550,ha='center',va='center')
    VL['ax'].legend(fontsize=8)

    # dt_object = datetime.fromtimestamp(ts_list[t]) # defaults to local time. Cool.
    dye_0lag['ax'].text(0.98,0.05,dt_list[t].astimezone(pdt).strftime("%m/%d/%Y,\n%H:%M PDT"),transform=dye_0lag['ax'].transAxes,ha='right',zorder=55,fontsize=8,color='k',fontweight='bold',bbox=dict(facecolor='white'))


    fig.subplots_adjust(right=0.98,left=0.1,bottom=0.15,top = 0.95,wspace=0.5)
    plt.savefig(home+f'plots/Oct 2024 samples/compare/station_profiles_{t:0>2}.png')
    # plt.show(block=False)
    # plt.pause(0.1)
    plt.close()


