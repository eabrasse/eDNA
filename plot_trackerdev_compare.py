import matplotlib.pyplot as plt
import string
import numpy as np
import netCDF4 as nc
from matplotlib.gridspec import GridSpec
import efun

plt.close('all')
tab10 = plt.get_cmap('tab10',10)

# load in data
home = '/Users/elizabethbrasseale/Projects/eDNA/'

data_fn = home+'data/trackercompare/release_2023.02.03.nc'
data_dev_fn = home+'data/trackercompare/release_2023.02.03_dev.nc'
data_dev2_fn = home+'data/trackercompare/release_2023.02.03_dev2.nc'
ds = nc.Dataset(data_fn)
ds_dev = nc.Dataset(data_dev_fn)
ds_dev2 = nc.Dataset(data_dev2_fn)


lon = {}
lon['xlabel'] = 'longitude (deg)'
lon['vname'] = 'lon'

lat = {}
lat['xlabel'] = 'latitude (deg)'
lat['vname'] = 'lat'

z = {}
z['xlabel'] = 'z (m)'
z['vname'] = 'z'

var_list = [lon,lat,z]
nvar = len(var_list)


times = [0,16]
ntimes = len(times)

fw,fh = efun.gen_plot_props()
fig = plt.figure(figsize=(ntimes*fw,2*fh))
gs = GridSpec(nvar,ntimes)

count = 0
for var in var_list:
    var['min'] = np.min([ds[var['vname']][:].min(),ds_dev[var['vname']][:].min(),ds_dev2[var['vname']][:].min()])
    var['max'] = np.max([ds[var['vname']][:].max(),ds_dev[var['vname']][:].max(),ds_dev2[var['vname']][:].max()])
    
    tcount = 0
    for time in times:
        ax = fig.add_subplot(gs[count,tcount])
        ax.hist(ds[var['vname']][time,:],bins=20,range=(var['min'],var['max']),facecolor=tab10(0),edgecolor=tab10(0),alpha=0.2)
        ax.hist(ds_dev[var['vname']][time,:],bins=20,range=(var['min'],var['max']),facecolor=tab10(1),edgecolor=tab10(1),alpha=0.2)
        ax.hist(ds_dev2[var['vname']][time,:],bins=20,range=(var['min'],var['max']),facecolor=tab10(2),edgecolor=tab10(2),alpha=0.2)
    
        ax.set_ylabel('count')
        ax.set_xlabel(var['xlabel'])
        
        if count==0:
            ax.text(0.1,0.9,f't = {time:} hours',transform=ax.transAxes)
            if tcount==0:
                ax.text(0.1,0.8,'tracker',transform=ax.transAxes,color=tab10(0))
                ax.text(0.1,0.7,'tracker_dev',transform=ax.transAxes,color=tab10(1))
                ax.text(0.1,0.6,'tracker_dev2',transform=ax.transAxes,color=tab10(2))
        
        tcount+=1
    
    count+=1

fig.subplots_adjust(top=0.98,left=0.15,right=0.98,bottom=0.1,hspace=0.3,wspace=0.3)
plt.show(block=False)
plt.pause(0.1)

ds.close()
ds_dev.close()