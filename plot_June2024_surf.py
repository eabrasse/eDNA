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

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])
tab10 = plt.get_cmap('tab10',10)
atoz = string.ascii_lowercase

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'

data_fn = home+'data/June2024_3dhist_w_ages.p'
D = pickle.load(open(data_fn,'rb'))
for key in D.keys():
    locals()[key] = D[key]
lon0 = -122.733651
lat0 = 47.740037

# gather mask from a random other data file, for plotting
data_fn = home+'data/HC_surfdye_0605.p'
Dp = pickle.load(open(data_fn,'rb'))
maskr = Dp['maskr'][:]

# load in sampling locations
sample_info_fn = home+'data/June2024_sampling_info.csv'
df = pd.read_csv(sample_info_fn,sep=',',engine='python')
Ds = {}
filter_type = df.filtered_with.unique()



for ft in filter_type:
    Ds[ft] = {}
    # I wanted to do this cleverly, but I think I need to hardcode it
    Ds[ft]['marker'] = 'o'
    Ds[ft]['markersize'] = 12
    Ds[ft]['mfc'] = tab10(0)
    Ds[ft]['zorder'] = 51
    if ft=='Ascension':
        Ds[ft]['marker'] = 'v'
        Ds[ft]['markersize'] = 8
        Ds[ft]['mfc'] = 'None'
        Ds[ft]['zorder']=55
    
    gg = df[df['filtered_with']==ft]
    Ds[ft]['sampling_locations'] = gg.sampling_location.unique()
    Ds[ft]['xsloc'] = np.zeros((len(Ds[ft]['sampling_locations'])))
    Ds[ft]['ysloc'] = np.zeros((len(Ds[ft]['sampling_locations'])))
    count =0
    for sampling_location in Ds[ft]['sampling_locations']:
        ggg = gg[gg['sampling_location']==sampling_location]
        # all data entries at each sampling location should have the same lat/lon
        # so just read the first index
        lon = ggg['long'].values[0]
        lat = ggg['lat'].values[0]
        Ds[ft]['xsloc'][count],Ds[ft]['ysloc'][count] = efun.ll2xy(lon,lat,lon0,lat0)
        count+=1

for t in range(len(ts_list)):
    # t=0 #while testing
    fw,fh = efun.gen_plot_props()
    fig = plt.figure(figsize=(fw*2,fh*1.75))
    ax = fig.gca()
    # gs = GridSpec(1,2)
    pmap = particle_map[t,-1,:,:]
    pmap[pmap==0]=np.nan
    ax.pcolormesh(x_edges,y_edges,maskr[1:-1,1:-1],cmap=cmap_mask,shading='flat',zorder=50)
    p=ax.pcolormesh(x_edges,y_edges,pmap,shading='flat',cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(vmin=10**0,vmax=10**3))
    
    for ft in filter_type:
        ax.plot(Ds[ft]['xsloc'],Ds[ft]['ysloc'],linestyle='None',marker=Ds[ft]['marker'],markersize=Ds[ft]['markersize'],markerfacecolor=Ds[ft]['mfc'],mec='k',zorder=Ds[ft]['zorder'])
    plt.colorbar(p)

    ax.set_xlabel('dist from pen (m)')
    ax.set_ylabel('dist from pen (m)')
    ax.set_aspect(1)

    dt_object = datetime.fromtimestamp(ts_list[t]) # defaults to local time. Cool.
    #dt_object = datetime.utcfromtimestamp(ts_list[t]) # if you want UTC
    ax.text(0.9,0.1,dt_object.strftime("%m/%d/%Y, %H:%M:%S PDT"),transform=ax.transAxes,ha='right',zorder=55)

    fig.subplots_adjust(right=0.98,left=0.125,bottom=0.08,top = 0.98,wspace=0.3)
    plt.savefig(home+f'plots/June 2024 samples/surface_{t}.png')
    # plt.show(block=False)
    # plt.pause(0.1)
    plt.close()


