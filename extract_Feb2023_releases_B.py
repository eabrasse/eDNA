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




#find Feb sampling particle tracks
track_dir0 = home+'LO_output/tracks/'

f_list = os.listdir(track_dir0)
f_list.sort()
f_list = [x for x in f_list if (x[:11]=='hc_dolph_3d')&(x[-4]!='6')&(x[-10:-6]=='2023')] # because releases were in Jan & Feb

# D = {}
# count=0
# print('Preprocessing particle release experiments')
# for f in f_list:
#     D[f] = {}
#
#     track_dir = track_dir0+f
#
#     # build a keyname from the release filename
#     file_list = os.listdir(track_dir)
#     file_list = [x for x in file_list if x[:3]=='rel']
#     rel_fn = file_list[0]
#
#     filefn = track_dir+'/'+rel_fn
#     D[f]['filefn']=filefn
#     ds = nc.Dataset(filefn)
#
#     ot = ds['ot'][:].data
#     T0 = datetime(1970,1,1,tzinfo=pytz.utc)+timedelta(seconds=ot[0]) #datetime of release
#     T1 = datetime(1970,1,1,tzinfo=pytz.utc)+timedelta(seconds=ot[-1])
#     D[f]['times'] = [T0,T1]
#
#     if count==0:
#         Tmin = T0
#         Tmax = T1
#         gridfn = track_dir + '/grid.nc'
#         dsg = nc.Dataset(gridfn)
#         lonr = dsg['lon_rho'][:]
#         latr = dsg['lat_rho'][:]
#         dsg.close()
#
#     else:
#         Tmin = min([Tmin,T0])
#         Tmax = max([Tmax,T1])
#
#     count+=1
#
# print('Preprocessing done!')
#set reference time
tz = pytz.utc
pdt = pytz.timezone('America/Vancouver')

# might... do this different
# dt_list0 = pd.date_range(start=Tmin,end = Tmax, freq="60min").to_pydatetime().tolist()
sample_dt = datetime(2023,2,3,13,0,tzinfo=pdt).astimezone(tz)

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

moor1 = {}
moor1['lon'] = -122.72845
moor1['lat'] = 47.747933
moor1['label'] = 'marginal pier'

moor2 = {}
moor2['lon'] = -122.729167
moor2['lat'] = 47.745633
moor2['label'] = 'between MP and bridge'


moor3 = {}
moor3['lon'] = -122.728933
moor3['lat'] = 47.74325
moor3['label'] = 'near bridge'

moor4 = {}
moor4['lon'] = -122.729783
moor4['lat'] = 47.742667
moor4['label'] = 'outside pool'

moor_list = [moor1,moor2,moor3,moor4]
nmoor = len(moor_list)


zref = -2

radius_list = np.arange(15,151,15)
nrad = len(radius_list)

for moor in moor_list:
    moor['count'] = np.zeros((nrad))

particle_count = np.zeros((nmoor,nrad))
particle_age_lists = [[[] for m in range(nmoor)] for r in range(nrad)]

for f in f_list:
    
    track_dir = track_dir0+f
    print(f'Working on {f}')


    # build a keyname from the release filename
    file_list = os.listdir(track_dir)
    file_list = [x for x in file_list if x[:3]=='rel']
    rel_fn = file_list[0]

    filefn = track_dir+'/'+rel_fn
    ds = nc.Dataset(filefn)
    ot = ds['ot'][:].data
    
    dt_list = [tz.localize(datetime(1970,1,1)+timedelta(seconds=ott)) for ott in ot]
    t = argnearest(dt_list,sample_dt)
    delta_T = ot[t]-ot[0]
    
    for m in range(nmoor):

        moor = moor_list[m]
        
        zmask = ds['z'][t,:]>(ds['zeta'][t,:]-2)
        
        if np.sum(zmask)>0: #(if there are any!)
        
            xp,yp = efun.ll2xy(ds['lon'][t,zmask],ds['lat'][t,zmask],moor['lon'],moor['lat']) # calculate x,y dist from mooring
            
            rp = np.sqrt(xp**2+yp**2) #calculate radial distance from mooring
            
            for r in range(nrad): #count how many are within different radiuses
                
                moor['count'][r] += np.sum(rp<radius_list[r])
                
                age_list = [delta_T]*int(np.sum(rp<radius_list[r]))
                
                for age in age_list:
                    particle_age_lists[r][m].append(age)

            
    ds.close()

print('Done extracting!')
print('Calculating mean age...')
particle_age_bins = np.zeros((nmoor,nrad))

for m in range(nmoor):
    for r in range(nrad):
        particle_age_bins[m,r] = np.mean(particle_age_lists[r][m])
print('Done calculating mean age!')

vnames = ['moor_list','particle_age_bins','radius_list','sample_dt','zref','particle_age_lists']
D={}
for vname in vnames:
    D[vname] = locals()[vname]
outfn = home+'LO_data/eDNA/Feb2023_moor_radius_counts.p'
pickle.dump(D,open(outfn, 'wb'))
print('saved to {}'.format(outfn))
