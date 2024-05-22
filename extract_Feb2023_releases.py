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
f_list = [x for x in f_list if (x[:11]=='hc_dolph_3d')&(x[-4]!='6')] # because releases were in Jan & Feb

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
        gridfn = track_dir + '/grid.nc'
        dsg = nc.Dataset(gridfn)
        lonr = dsg['lon_rho'][:]
        latr = dsg['lat_rho'][:]
        dsg.close()
        
    else:
        Tmin = min([Tmin,T0])
        Tmax = max([Tmax,T1])
    
    count+=1

print('Preprocessing done!')
#set reference time
Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

# might... do this different
dt_list0 = pd.date_range(start=Tmin,end = Tmax, freq="60min").to_pydatetime().tolist()

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773


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

xbins = 99
ybins = 101
bin_lon_edges=np.linspace(aa[0], aa[1],xbins+1)
bin_lat_edges=np.linspace(aa[2], aa[3],ybins+1)
bin_x_edges=np.linspace(aax[0], aax[1],xbins+1)
bin_y_edges=np.linspace(aay[0], aay[1],ybins+1)


nt = len(dt_list0)
hist_particles_bin = np.zeros((nt,ybins,xbins))
particle_age_lists = [[[[] for i in range(xbins)] for j in range(ybins)] for t in range(nt)]
total_released_particles = np.zeros((nt))
deep_particles = np.zeros((nt))
outfield_particles = np.zeros((nt))

count =0
flag=0
for dt in dt_list0:
    
    print('Extracting particle positions from all releases at timestep {}'.format(dt))
    
    
    for f in f_list:

        [T0,T1] = D[f]['times']

        if (dt>T0)&(dt<T1):

            
            # open actual particle tracks
            ds = nc.Dataset(D[f]['filefn'])
            
            if flag==0:
                particle_rel=np.shape(ds['lon'])[1]
                flag+=1
                
            # add active particles
            total_released_particles[count] += particle_rel
            
            # find specific time index
            ot = ds['ot'][:].data
            
            dt_list = [utc.localize(datetime(1970,1,1)+timedelta(seconds=ott)) for ott in ot]
            t = argnearest(dt_list,dt)
            delta_T = ot[t]-ot[0] # calculate time since release
            
            #mask particles deeper than 2m below surface
            zmask = ds['z'][t,:]>(ds['zeta'][t,:]-2)
            
            # count how many won't be in the histogram because they're depth masked
            deep_particles[count] += np.sum(zmask)
            
            # now look at the surface particles...
            if np.sum(zmask)>0: #(if there are any!)
                #grab x/y locations for surface particles at time stamp t 
                xp,yp = efun.ll2xy(ds['lon'][t,zmask],ds['lat'][t,zmask],lon0,lat0)
                # count how many total are in the surface layer (should be totalreleased - zmasked)
                all_surface_particles[count] += len(xp.flatten())
                # bin them by location
                hist = np.histogram2d(xp,yp,bins=[bin_x_edges,bin_y_edges])
                
                #some will be outside the histogram domain. Estimate how many these are by subtracting the binned ones from all surface ones
                outfield_particles[count] += len(xp.flatten())-np.sum(hist[0])
                
                #save the histogram
                hist_particles_bin[count,:] += hist[0]
                
                for i in range(xbins):
                    for j in range(ybins):
                        age_list = [deltaT]*int(hist[0][j,i])
                        particle_age_lists[count][j][i].append(age_list)
                
            ds.close()
    count+=1
print('Done extracting!')
print('Calculating mean age...')
particle_age_bins = np.zeros((nt,ybins,xbins))
for t in range(nt):
    for i in range(xbins):
        for j in range(ybins):
            particle_age_bins[t,j,i] = np.mean(particle_age_lists[t][j][i])
print('Done calculating mean age!')

vnames = ['hist_particles_bin','outfield_particles','deep_particles','particle_age_bins','bin_x_edges','bin_y_edges','dt_list0','zref']
D={}
for vname in vnames:
    D[vname] = locals()[vname]
outfn = home+'LO_data/eDNA/Feb2023_hist_counts.p'
pickle.dump(moor_dict,open(outfn, 'wb'))
print('saved to {}'.format(outfn))
