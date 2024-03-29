import os
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from PIL import Image
import numpy as np
from datetime import datetime, timedelta
import pytz
import netCDF4 as nc
import pandas as pd
import pickle
import efun

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)

plt.close('all')
home = '/data2/pmr4/eab32/'

#find June sampling particle tracks
track_dir0 = home+'LO_output/tracks/'

f_list = os.listdir(track_dir0)
f_list.sort()
f_list = [x for x in f_list if (x[:11]=='hc_dolph_3d')&(x[-4]=='9')] # because releases were in Sep

scale_const = False
TV = False

D = {}
count=0
print('Preprocessing particle release experiments')
for f in f_list:
    D[f] = {}

    track_dir = track_dir0+f

    # build a keyname from the release filename
    file_list = os.listdir(track_dir)
    file_list = [x for x in file_list if x[:3]=='rel']
    rel_fn = file_list[0]

    filefn = track_dir+'/'+rel_fn
    D[f]['filefn']=filefn
    ds = nc.Dataset(filefn)

    ot = ds['ot'][:].data
    T0 = datetime(1970,1,1,tzinfo=pytz.utc)+timedelta(seconds=ot[0]) #datetime of release
    T1 = datetime(1970,1,1,tzinfo=pytz.utc)+timedelta(seconds=ot[-1])
    D[f]['times'] = [T0,T1]
    
    if count==0:
        Tmin = T0
        Tmax = T1
        if scale_const:
            DNA_mean = np.mean(df.PB_quantity_mean)
        gridfn = track_dir + '/grid.nc'
        dsg = nc.Dataset(gridfn)
        lonr = dsg['lon_rho'][:]
        latr = dsg['lat_rho'][:]
        dsg.close()
        
    else:
        Tmin = min([Tmin,T0])
        Tmax = max([Tmax,T1])
    if TV:
        D[f]['C0'] = df[df.date==T0].PB_quantity_mean.values[0] 
    
    count+=1

print('Preprocessing done!')
#set reference time
Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

# first_sample_utc = Pacific.localize(datetime(2023,1,31,0,0)).astimezone(utc)
# first_release = first_sample_utc-timedelta(days=2)
# last_sample_utc = Pacific.localize(datetime(2023,2,2,0,0)).astimezone(utc)
# might... do this different
dt_list0 = pd.date_range(start=Tmin,end = Tmax, freq="60min").to_pydatetime().tolist()

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

moor0 = {}
moor0['lon'] = lon0
moor0['lat'] = lat0
moor0['label'] = 'Dolphin pen'

moor1 = {}
moor1['lon'] = -122.729738
moor1['lat'] = 47.74284
moor1['label'] = 'Inside dolphin pen - 3 dolphins'

moor2 = {}
moor2['lon'] = -122.729626
moor2['lat'] = 47.742852
moor2['label'] = 'outside dolphin pen'

moor3 = {}
moor3['lon'] = -122.728678
moor3['lat'] = 47.748029
moor3['label'] = 'marginal pier'

moor4 = {}
moor4['lon'] = -122.73002
moor4['lat'] = 47.742233
moor4['label'] = 'near sea lions'

moor5 = {}
moor5['lon'] = -122.729785
moor5['lat'] = 47.742766
moor5['label'] = 'Inside dolphin pen - 1 dolphin'

moor_list = [moor0,moor1,moor2,moor3,moor4,moor5]


# MAP
xgrid,ygrid = efun.ll2xy(lonr,latr,lon0,lat0)
# set domain limits
pad = .01
# plot region around delta pier, using release point as indicator
aa = [lon0-pad, lon0+pad,
   lat0-pad, lat0+pad]

aax,aay =  efun.ll2xy(np.array(aa[:2]),np.array(aa[2:]),lon0,lat0)
aaxy = [aax[0],aax[1],aay[0],aay[1]]

#identify grid edge limits for making mask      
AA = [lonr[0,0], lonr[0,-1],
        latr[0,0], latr[-1,0]]
zref = -2

nbins = 100
bin_lon_edges=np.linspace(aa[0], aa[1],nbins+1)
bin_lat_edges=np.linspace(aa[2], aa[3],nbins+1)
bin_x_edges=np.linspace(aax[0], aax[1],nbins+1)
bin_y_edges=np.linspace(aay[0], aay[1],nbins+1)
# xx, yy = np.meshgrid(bin_x_edges[:-1]+0.5*(bin_x_edges[1]-bin_x_edges[0]),bin_y_edges[:-1]+0.5*(bin_y_edges[1]-bin_y_edges[0]))
dxbin = bin_x_edges[1]-bin_x_edges[0]
dybin = bin_y_edges[1]-bin_y_edges[0] 
dzbin = np.abs(zref)
vol_bin = dxbin*dybin*dzbin

# find which bin
for moor in moor_list:
    #first bin they're greater than
    moor['lon_bin'] = np.argwhere(moor['lon']<bin_lon_edges)[0][0]
    moor['lat_bin'] = np.argwhere(moor['lat']<bin_lat_edges)[0][0]

nt = len(dt_list0)
const_DNA_bin = np.zeros((nt,nbins,nbins))
active_particles = np.zeros((nt))

k_decay = 0.02/3600 #units: data 0.02 1/hr, multiply by hr/sec to get 1/sec

count =0
flag=0
for dt in dt_list0:
    
    print('Extracting particle positions from all releases at timestep {}'.format(dt))
    
    
    
    #FIGURE OUT HOW TO INCORPORATE TIME VARYING WITHOUT LOOKING IT UP OVER AND OVER
    for f in f_list:

        [T0,T1] = D[f]['times']

        # ot = ds['ot'][:].data
#
#         ot0 = utc.localize(datetime(1970,1,1)+timedelta(seconds=ot[0]))
#         ot1 = utc.localize(datetime(1970,1,1)+timedelta(seconds=ot[-1]))

        if (dt>T0)&(dt<T1):
            ds = nc.Dataset(D[f]['filefn'])
            
            if flag==0:
            #     xp0,yp0 = efun.ll2xy(ds['lon'][0,:],ds['lat'][0,:],lon0,lat0)
            #     dxrel = np.abs(xp0.max()-xp0.min())
            #     dyrel = np.abs(yp0.max()-yp0.min())
            #     dzrel = 1 #all released at surface, soooo....
            #     vol_rel = dxrel*dyrel*dzrel
            #
                particle_rel=np.shape(ds['lon'])[1]
            #     particle_conc_rel = particle_rel/vol_rel
            #
                flag+=1
                
            ot = ds['ot'][:].data
            
            dt_list = [utc.localize(datetime(1970,1,1)+timedelta(seconds=ott)) for ott in ot]
            t = argnearest(dt_list,dt)
            delta_T = ot[t]-ot[0]
            decay = np.exp(-k_decay*delta_T)
            zmask = ds['z'][t,:]>(ds['zeta'][t,:]-2)
            
            active_particles[count] += particle_rel
            
            if np.sum(zmask)>0:
                xp,yp = efun.ll2xy(ds['lon'][t,zmask],ds['lat'][t,zmask],lon0,lat0)
                # ax.scatter(ds['lon'][t,zmask],ds['lat'][t,zmask],c='w',s=1,alpha=0.05)
                # hist = np.histogram2d(ds['lon'][t,zmask],ds['lat'][t,zmask],bins=[bin_lon_edges,bin_lat_edges])
                hist = np.histogram2d(xp,yp,bins=[bin_x_edges,bin_y_edges])
                
                particle_conc_bin = hist[0]/vol_bin
                # DNA_conc_bin = DNA_conc_rel * particle_conc_bin/particle_conc_rel
                
                #note: not tranposing because I'm not plotting in 2d! I want to be able to follow the
                # indexing in the manual. Histogram should be indexed nx, ny
                # const_particle_bin[count,:] += decay*hist[0]
                # TV_particle_bin[count,:] += D[f]['C0']*decay*hist[0]
                if scale_const:
                    const_DNA_bin[count,:] += DNA_mean*decay*particle_conc_bin/particle_conc_rel
                else:
                    const_DNA_bin[count,:] += decay*particle_conc_bin#/particle_conc_rel

                
                
            ds.close()
    count+=1

print('Done extracting!')
print('Extracting moorings')
for moor in moor_list:
    # moor['const_particle_bin'] = const_particle_bin[:,moor['lon_bin'],moor['lat_bin']]
    # moor['TV_particle_bin'] = TV_particle_bin[:,moor['lon_bin'],moor['lat_bin']]
    moor['const_DNA_bin'] = const_DNA_bin[:,moor['lon_bin'],moor['lat_bin']]

print('Done!')

moor_dict = {'moor_list':moor_list,'const_DNA_bin':const_DNA_bin,'dt_list':dt_list0,'active_particles':active_particles,'bin_x_edges':bin_x_edges,'bin_y_edges':bin_y_edges}


outfn = home+'LO_data/eDNA/Sep2022_DNA_moorings_extended_const_only_noscale.p'
pickle.dump(moor_dict,open(outfn, 'wb'))
print('saved to {}'.format(outfn))
