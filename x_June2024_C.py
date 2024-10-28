import os
import sys
import numpy as np
from datetime import datetime, timedelta
import pytz
import netCDF4 as nc
import pandas as pd
import pickle
import efun
from lo_tools import zrfun

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

#set dt list based on samples
Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

home = '/data2/pmr4/eab32/'

sample_info_fn = home+'LO_data/eDNA/June2024_sampling_info.csv'
df = pd.read_csv(sample_info_fn,sep=',',engine='python')
old_columns = [col for col in df.columns]

df['dt0'] = pd.to_datetime(df['date_filtered']+' '+df['time_filtered_t0']).dt.tz_localize(Pacific)
df['ts0'] = df.dt0.values.astype(np.int64) // 10 ** 9
# t00 = np.nanmin(df.ts0.values[:])
# t11 = np.nanmax(df.ts0.values[:])
# df['ts_norm'] = (df.ts0.values[:]-t00)/(t11-t00) # scale from 0 to 1
xsloc,ysloc = efun.ll2xy(df['long'].values[:],df['lat'].values[:],lon0,lat0)
df['xsloc'] = xsloc
df['ysloc'] = ysloc

#find June sampling particle tracks
track_dir0 = home+'LO_output/tracks2/hc11_v01_uu0k/'

f_list = os.listdir(track_dir0)
f_list.sort()
f_list = [x for x in f_list if (x[:11]=='hc_dolph_3d')&(x[17:24]=='2024.06')]
# note current these are the only releases I've done with the new tracker code, so no need to subselect



first_sample_utc = Pacific.localize(datetime(2024,6,13,10,0)).astimezone(utc)
last_sample_utc = Pacific.localize(datetime(2024,6,13,17,0)).astimezone(utc)
dt_list0 = pd.date_range(start=first_sample_utc,end = last_sample_utc, freq="60min").to_pydatetime().tolist()
ts_list = [datetime.timestamp(dt) for dt in dt_list0]
dt = ts_list[1]-ts_list[0]
nt = len(ts_list)


df['t0i'] = np.zeros(len(df.index)).astype(int)
df['ts0i'] = np.zeros(len(df.index))
df['t1i'] = np.zeros(len(df.index)).astype(int)
df['ts1i'] = np.zeros(len(df.index))
for ind in df.index:
    df['t0i'][ind] = int(np.argwhere(ts_list<df['ts0'][ind])[-1][0])
    df['ts0i'][ind] = ts_list[df['t0i'][ind]]
    if np.sum(ts_list>df['ts0'][ind])>0:
        df['t1i'][ind] = int(np.argwhere(ts_list>df['ts0'][ind])[0][0])
        df['ts1i'][ind] = ts_list[df['t1i'][ind]]
    else:
        df['t1i'][ind] = np.nan
        df['ts1i'][ind] = np.nan
# df['t0t1'] = [[t0i[i],t1i[i]] for i in range(len(t0i))]

his_dir = '/data1/jxiong/LO_roms/hc11_v01_uu0k/'
his_fn_list = [his_dir + 'f2024.06.13/ocean_his_0018.nc',
                his_dir + 'f2024.06.13/ocean_his_0019.nc',
                his_dir + 'f2024.06.13/ocean_his_0020.nc',
                his_dir + 'f2024.06.13/ocean_his_0021.nc',
                his_dir + 'f2024.06.13/ocean_his_0022.nc',
                his_dir + 'f2024.06.13/ocean_his_0023.nc',
                his_dir + 'f2024.06.13/ocean_his_0024.nc',
                his_dir + 'f2024.06.13/ocean_his_0025.nc']
his_fn_list.sort()

tt=0
for fn in his_fn_list:
    ds = nc.Dataset(fn)
    if fn==his_fn_list[0]:
        nz,ny,nx = ds['salt'][0,:,:,:].shape
        z_edges = np.zeros((len(his_fn_list),nz,ny,nx))
        h = ds['h'][:]
        G,S,T = zrfun.get_basic_info(fn)
        lonr = G['lon_rho'][:]
        latr = G['lat_rho'][:]
        xr,yr = efun.ll2xy(lonr,latr,lon0,lat0)
        df['xli'] = df.apply(lambda row: np.argmin(np.abs(xr[0,:]-row.xsloc)), axis=1)
        df['yli'] = df.apply(lambda row: np.argmin(np.abs(yr[:,0]-row.ysloc)), axis=1)
    zeta = ds['zeta'][0,:]
    z_edges [tt,:] = zrfun.get_z(h, zeta, S, only_rho=True)

    ds.close()
    tt+=1


pz0 = np.zeros((len(df.index),nz-1))
pz1 = np.zeros((len(df.index),nz-1))
# particle_age_lists = [[[[[] for x in range(nx-2)] for y in range(ny-2)] for z in range(nz-2)] for t in range(nt)]
rad = 100
depth = 5

# what I want to do is assign each station the before & after time step,
# then at each time step identify which stations are "active" at those time steps
# and loop through active stations & count up particles
# maybe I should save particles at t0 and t1 during the loops
# then interpolate between t0 and t1 afterwards


count=1

# for f in f_list:

f=f_list[0] #testing

print(f'working on file {count} of {len(f_list)}')

track_dir = track_dir0+f

# build a keyname from the release filename
file_list = os.listdir(track_dir)
file_list = [x for x in file_list if x[:3]=='rel']
rel_fn = file_list[0]

filefn = track_dir+'/'+rel_fn
ds = nc.Dataset(filefn)

ot = ds['ot'][:].data
ts_list_p = []
for tt in ot:
    ts_list_p.append(datetime.timestamp(datetime(1970,1,1,tzinfo=pytz.utc)+timedelta(seconds=tt)))

for t in range(nt):
    
    print(f'Time step {t}...')
    dt_list = [np.abs(ts_p-ts_list[t]) for ts_p in ts_list_p]
    pt = np.argmin(dt_list)
    xp,yp = efun.ll2xy(ds['lon'][pt,:],ds['lat'][pt,:],lon0,lat0)
    # pt = np.argmin(np.abs(ts_list_p-ts_list[t]))
    # delta_T = ts_list_p[pt]-ts_list_p[0]
    
    # count = np.sum((np.abs(xp-df.xsloc)<rad)*(np.abs(yp-df.ysloc)<rad)*(np.abs(ds['z'][pt,:]-df.depth_m)<depth))
    ind0 = np.argwhere(df.ts0i==ts_list[t])[0][:]
    for ind in ind0:
        rpm = np.sqrt((xp-df.xsloc[ind])**2+(yp-df.ysloc[ind])**2)<100
        count,edges = np.histogram(ds['z'][pt,rpm],z_edges[t,:,df.yli[ind],df.xli[ind]])
        pz0[ind,:] += count[:]
        
    ind1 = np.argwhere(df.ts1i==ts_list[t])[0][:]
    for ind in ind1:
        rpm = np.sqrt((xp-df.xsloc[ind])**2+(yp-df.ysloc[ind])**2)<100
        count,edges = np.histogram(ds['z'][pt,rpm],z_edges[t,:,df.yli[ind],df.xli[ind]])
        pz1[ind,:] += count[:]
    

ds.close()
count+=1

particle_profiles = np.zeros((len(df.index),nz-1))
for ind in df.index:
    # if np.isnan(df[ind].ts)
    particle_profiles[ind,:] = pz0[ind,:]+(df.ts0[ind]-df.ts0i[ind])*(pz1[ind,:]-pz0[ind,:])/dt
# df['ps'] = df.apply(lambda row: row.p0 + (row.ts0-row.ts0i)*(row.p1-row.p0)/dt, axis=1)

D = {}
var_list = ['particle_profiles','z_edges','ts_list']
for var in var_list:
    D[var] = locals()[var]


outfn = home+'LO_data/eDNA/June2024_sample_profiles.p'
pickle.dump(D,open(outfn, 'wb'))
print('saved to {}'.format(outfn))
