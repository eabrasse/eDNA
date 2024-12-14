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
xx,yy = np.meshgrid(x_center,y_center)
dx0 = np.diff(x_edges) #not uniform
dy0 = np.diff(y_edges) # not uniform
dx,dy = np.meshgrid(dx0,dy0)
dz = z_edges[1]-z_edges[0] #uniform
nz = len(z_center)
nt = len(ts_list)
dv0 = dx*dy*dz #2D array
dv=np.tile(dv0,(nt,nz,1,1)) # 4d array
# particle_conc_map = particle_map/dv

dt = ts_list[1]-ts_list[0]
# particle_map = particle_map/dv

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

fw,fh = efun.gen_plot_props()
fig = plt.figure(figsize=(fw*2,fh*1.25))
ax = fig.gca()
ps_list_20m = []
ps_list_100m = []
label20m = 'xy = HC grid, z = 10m'
label100m = 'xy = 100 m rad, z = 10m'
for ind in df.index:
    xi = np.argmin(np.abs(x_center-df['xsloc'][ind]))
    yi = np.argmin(np.abs(y_center-df['ysloc'][ind]))
    zi = np.argmin(np.abs(z_center-df['depth_m'][ind]))
    
    # now estimate 100m
    dr = np.sqrt((xx-df['xsloc'][ind])**2+(yy-df['ysloc'][ind])**2)
    
    t0i = np.argwhere(ts_list<df['ts0'][ind])[-1][0] #should be max
    p0tz = particle_map[t0i,zi,:]
    p0 = p0tz[yi,xi]/dv0[yi,xi]
    p0_100 = np.sum(p0tz[dr<100])/np.sum(dv0[dr<100])
    if np.sum(ts_list>df['ts0'][ind])>0:
        
        t1i = np.argwhere(ts_list>df['ts0'][ind])[0][0] # should be min
        p1tz = particle_map[t1i,zi,:]
        
        p1 = p1tz[yi,xi]/dv0[yi,xi]
        ps = p0 + (df['ts0'][ind]-ts_list[t0i])*(p1-p0)/dt
        
        p1_100 = np.sum(p1tz[dr<100])/np.sum(dv0[dr<100])
        ps_100 = p0_100 + (df['ts0'][ind]-ts_list[t0i])*(p1_100-p0_100)/dt
        
    else:
        t1i=t0i
        
        p1=p0
        ps = p0
        
        p1=p0_100
        ps_100 = p0_100
        
    ps_list_20m.append(ps)
    ps_list_100m.append(ps_100)
    
    ax.plot([ts_list[t0i],ts_list[t1i]],[p0,p1],color=tab10(0),marker='none',linestyle='solid',lw=1.5)
    ax.plot([ts_list[t0i],ts_list[t1i]],[p0_100,p1_100],color=tab10(1),marker='none',linestyle='solid',lw=1.5)
    ax.axvline(x=df['ts0'][ind],color='gray',linestyle='dashed',lw=1.5)
    
    if ind>0:
        label20m='_nolegend_'
        label100m='_nolegend_'
    ax.plot(df['ts0'][ind],ps,marker='o',markersize=8,mec='None',mfc=tab10(0),alpha=0.5,label=label20m)
    ax.plot(df['ts0'][ind],ps_100,marker='o',markersize=8,mec='None',mfc=tab10(1),alpha=0.5,label=label100m)


ax.set_xlabel('Sample timestamp (seconds)')
ax.set_ylabel('Particle concentration')
# ax.set_yscale('log')
ax.legend()

plabel = 'Model (particles m-3)'
df['Model (particles m-3)'] = ps_list_100m

old_columns.append(plabel)
df = df[old_columns]
new_sample_info_fn = home+'data/June2024_sampling_info_w_pred.csv'
df.to_csv(new_sample_info_fn)


# fig.subplots_adjust(right=0.98,left=0.05,bottom=0.11,top = 0.98,wspace=0.3)
plt.savefig(home+f'plots/June 2024 samples/predictions.png')
# plt.show(block=False)
# plt.pause(0.1)
plt.close()




