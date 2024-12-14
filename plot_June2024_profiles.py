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

data_fn = home+'data/June2024_3dhist_zw_no_ages.p'
D = pickle.load(open(data_fn,'rb'))
for key in ['x_edges','y_edges','z_edges','ts_list','particle_map']:
    locals()[key] = D[key]

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773
x_center = 0.5*(x_edges[1:]+x_edges[:-1])
y_center = 0.5*(y_edges[1:]+y_edges[:-1])
z_center = 0.5*(z_edges[:,1:,:,:]+z_edges[:,:-1,:,:])
z_center = 0.5*(z_center[:,:,:,1:]+z_center[:,:,:,:-1])
z_center = 0.5*(z_center[:,:,1:,:]+z_center[:,:,:-1,:])
xx,yy = np.meshgrid(x_center,y_center)
dx = np.diff(x_edges) #not uniform
dy = np.diff(y_edges) # not uniform
# dx,dy = np.meshgrid(dx0,dy0)
yl0 = z_center.min()

# dz = z_edges[1]-z_edges[0] #uniform
dz = np.diff(z_edges,axis=1)
dz = 0.5*(dz[:,:,:,1:]+dz[:,:,:,:-1])
dz = 0.5*(dz[:,:,1:,:]+dz[:,:,:-1,:])
nt,nz,ny,nx = np.shape(dz)
dx = np.tile(np.reshape(dx,[1,1,1,nx]),[nt,nz,ny,1])
dy = np.tile(np.reshape(dy,[1,1,ny,1]),[nt,nz,1,nx])
dv = dx*dy*dz #4d array

dt = ts_list[1]-ts_list[0]

zw_fn = home+'data/June2024_sample_profiles.p'
Dz = pickle.load(open(zw_fn,'rb'))
particle_profiles = Dz['particle_profiles'][:]
particle_profiles = 0.5*(particle_profiles[:,1:]+particle_profiles[:,:-1])


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
old_columns = [col for col in df.columns]

df['dt0'] = pd.to_datetime(df['date_filtered']+' '+df['time_filtered_t0']).dt.tz_localize(pdt)
df['ts0'] = df.dt0.values.astype(np.int64) // 10 ** 9
# t00 = np.nanmin(df.ts0.values[:])
# t11 = np.nanmax(df.ts0.values[:])
# df['ts_norm'] = (df.ts0.values[:]-t00)/(t11-t00) # scale from 0 to 1
xsloc,ysloc = efun.ll2xy(df['long'].values[:],df['lat'].values[:],lon0,lat0)
df['xsloc'] = xsloc
df['ysloc'] = ysloc
xmin = np.nanmin(xsloc)
xmax = np.nanmax(xsloc)
ymin = np.nanmin(ysloc)
ymax = np.nanmax(ysloc)
axis_fudge=500
a_ll = [xmin-axis_fudge,xmax+axis_fudge,ymin-axis_fudge,ymax+axis_fudge]


particle_conc = np.zeros(np.shape(particle_profiles))
for ind in range(np.shape(particle_profiles)[0]):
    ti = np.argmin([np.abs(ts-df.ts0[ind]) for ts in ts_list])
    xi = np.argmin(np.abs(x_center-df.xsloc[ind]))
    yi = np.argmin(np.abs(y_center-df.ysloc[ind]))
    zi = np.argmin(np.abs(z_center[ti,:,yi,xi]+df.depth_m[ind]))
    particle_conc[ind,:] = particle_profiles[ind,:]/dv[ti,:,yi,xi]
xl1 = particle_conc.max()
# Ds = {}
# filter_type = df.filtered_with.unique()
# for ft in filter_type:
#     Ds[ft] = {}
#     # I wanted to do this cleverly, but I think I need to hardcode it
#     Ds[ft]['marker'] = 'o'
#     Ds[ft]['markersize'] = 12
#     Ds[ft]['mfc'] = 'None'
#     Ds[ft]['zorder'] = 75
#     if ft=='Ascension':
#         Ds[ft]['marker'] = 'v'
#         Ds[ft]['markersize'] = 8
#         Ds[ft]['mfc'] = 'None'
#         Ds[ft]['zorder']=80


# ndeep = len(df[df.filtered_with=='Ascension'].sampling_location.unique()) # number of deep stations

# axis_fudge=500



for ind in df.index:
    if np.isnan(df.xsloc[ind]):
        continue
        
    # ind = 30 #testing
    fw,fh = efun.gen_plot_props()
    fig = plt.figure(figsize=(fw*1.3,fh*1.25))
    ax = fig.gca()

    ti = np.argmin([np.abs(ts-df.ts0[ind]) for ts in ts_list])
    xi = np.argmin(np.abs(x_center-df.xsloc[ind]))
    yi = np.argmin(np.abs(y_center-df.ysloc[ind]))
    zi = np.argmin(np.abs(z_center[ti,:,yi,xi]+df.depth_m[ind]))

    # plot profile of particle concentrations
    ax.plot(particle_conc[ind,:],z_center[ti,:,yi,xi],linestyle='solid',marker='none',color='k',lw=2)
    # plot zero line
    # ax.axvline(x=0,linestyle='dashed',lw=0.5,color='k')
    # plot sample depth & prediction
    ax.plot(particle_conc[ind,zi],z_center[ti,zi,yi,xi],linestyle='none',marker='o',mfc=tab10(6),mec='k',markersize=8)

    ax.set_xlabel('Particle concentration (m-3)')
    ax.set_ylabel('Depth (m)')
    # ax.set_xscale('log')
    # if ind==0:
    #     x0,x1 = ax.get_xlim()
    #     y0,y1 = ax.get_ylim()
    ax.set_xlim(-0.1,xl1+0.5)
    ax.set_ylim(yl0,5)
    ax.spines[['right', 'top']].set_visible(False)

    #add context map
    ax_in = inset_axes(ax, width="30%", height="30%", loc='upper right',bbox_transform=ax.transAxes,bbox_to_anchor=(-0.12,-0.06,1,1))
    
    zmask = z_center[ti,:,:,:]>-10
    pmap = np.sum(zmask*particle_map[ti,:,:,:],axis=0)/np.sum(zmask*dv[ti,:,:,:],axis=0)
    pmap[pmap==0]=np.nan
    
    # land mask
    ax_in.pcolormesh(xr,yr,maskr,cmap=cmap_mask,shading='nearest',zorder=45)
    # bathymetry
    # p=ax_in.pcolormesh(xr,yr,-h,cmap=cmo.cm.deep_r,shading='nearest',zorder=44)
    # particle field
    p=ax_in.pcolormesh(x_center,y_center,pmap,shading='nearest',cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(vmin=10**-3,vmax=3*10**-1))
    # bathy contours
    ax_in.contour(xr,yr,maskr,levels=[0.5],colors='k',linewidths=0.8,zorder=50,linestyles='-')
    ax_in.contour(xr,yr,h,levels=np.arange(0,150,10),colors='k',linewidths=0.4,zorder=50,linestyles='-',alpha=0.5)
    # sampling location
    ax_in.plot(df.xsloc[ind],df.ysloc[ind],linestyle='none',marker='o',mfc='k',mec='None',markersize=1,zorder=50)
    # area averaged over to make prediction
    circle = plt.Circle((df.xsloc[ind], df.ysloc[ind]), 100, color='k', fill=False,hatch='x',lw=0.8)
    ax_in.add_patch(circle)
    
    ax_in.set_xlabel('Dist (m)')
    ax_in.set_ylabel('Dist (m)')
    ax_in.set_aspect(1)
    ax_in.axis(a_ll)
    cbaxes = inset_axes(ax_in, width="4%", height="40%", loc='center right',bbox_transform=ax_in.transAxes,bbox_to_anchor=(0.,-0.1,1,1))
    cbaxes.tick_params(axis='both',which='both',labelsize=8,size=2)
    cb = fig.colorbar(p, cax=cbaxes,orientation='vertical')
    cb.set_label('particles m-3',fontsize=8)

    # add info
    refstr = 'lab ID: {}\nreferfence: {}\nsample time: {}\npred particle conc: {:.3} m-3'.format(df.lab_ID[ind],df.reference[ind],df.dt0[ind].strftime("%m/%d/%Y, %H:%M PDT"),particle_conc[ind,zi])
    ax.text(0.98,0.02,refstr,ha='right',va='bottom',transform=ax.transAxes,fontsize=10)

    # fig.subplots_adjust(right=0.98,left=0.05,bottom=0.11,top = 0.98,wspace=0.3)
    fig.tight_layout()
    plt.savefig(home+f'plots/June 2024 samples/sample_pmap_{ind:0>2}.png')
    # plt.show(block=False)
    # plt.pause(0.1)
    plt.close()




