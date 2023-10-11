import os
import sys
import numpy as np
from datetime import datetime, timedelta
import pytz
import netCDF4 as nc
import pandas as pd
import efun
import pickle
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmocean as cmo
from matplotlib.gridspec import GridSpec
import string


tab10 = plt.get_cmap('tab10',10)
atoz= string.ascii_lowercase

home = '/data2/pmr4/eab32/'

# get grid from a random file
data_fn = home+'LO_output/tracks/all_3d_hc_dolph_releases.p'

D = pickle.load(open(data_fn,'rb'))

# gather some fields, for convenience
lonp = D['metadata']['lon'][:]
latp = D['metadata']['lat'][:]
hh = D['metadata']['h'][:]
maskr = D['metadata']['mask'][:]

#get mean concentration 
data_fn = home+'LO_data/eDNA/ESP_Feb2023_hourly.csv'
df = pd.read_csv(data_fn,sep=',',engine='python')
DNA_mean = np.mean(df.PB_quantity_mean)

# load mooring results
moor_fn = home+'LO_data/eDNA/Feb2023_moorings.p'

moor_dict = pickle.load(open(moor_fn,'rb'))
const= {}
const['particle_bin'] = DNA_mean*moor_dict['const_particle_bin']
const['label'] = 'Constant'
TV = {}
TV['particle_bin'] = moor_dict['TV_particle_bin']
TV['label'] = 'Time-varying'
relstrat_list = [const,TV]
for relstrat in relstrat_list:
    relstrat['mean'] = np.nanmean(relstrat['particle_bin'],axis=0)
    relstrat['std'] = np.nanstd(relstrat['particle_bin'],axis=0)
    relstrat['var'] = relstrat['std']**2
    relstrat['cv'] = relstrat['std']/relstrat['mean']


#didn't include my stupid bins - here they are
# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773
pad = .01
aa = [lon0-pad, lon0+pad,
   lat0-pad, lat0+pad]
nbins = 100
bin_lon_edges=np.linspace(aa[0], aa[1],nbins+1)
bin_lat_edges=np.linspace(aa[2], aa[3],nbins+1)
xx, yy = np.meshgrid(bin_lon_edges[:-1]+0.5*(bin_lon_edges[1]-bin_lon_edges[0]),bin_lat_edges[:-1]+0.5*(bin_lat_edges[1]-bin_lat_edges[0]))

statfig_list = ['mean','var','cv']
norm_list = {'mean':matplotlib.colors.LogNorm(vmin=2e1,vmax=2e4),'var':matplotlib.colors.LogNorm(vmin=2e2,vmax=2e9),'cv':None}
colmap = cmo.cm.matter
for stat in statfig_list:
    
    plt.close('all')
    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(2,len(relstrat_list)+1)
    axlon = fig.add_subplot(gs[1,len(relstrat_list)])
    axlat = fig.add_subplot(gs[0,len(relstrat_list)])
    
    rscount=0
    for relstrat in relstrat_list:
        axmap = fig.add_subplot(gs[:,rscount])
        axmap.axis(aa)
        axmap.contour(lonp,latp,maskr,levels=[0.5],colors='k',linewidths=1,linestyles='solid')
    
            
        p=axmap.pcolormesh(xx,yy,relstrat[stat].T,cmap = colmap,norm= norm_list[stat]) 
        cbaxes = inset_axes(axmap, width="4%", height="40%", loc='center right',bbox_transform=axmap.transAxes,bbox_to_anchor=(-0.2,-0.2,1,1))
        cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
        cb.set_label('{} particle wt'.format(stat))
        
        axmap.text(0.1,0.9,'{}) {} {}'.format(atoz[rscount],relstrat['label'],stat),transform=axmap.transAxes)
        axmap.set_xlabel('Longitude')
        axmap.set_ylabel('Latitude')
        
        axlon.scatter(xx[0,:],np.mean(relstrat[stat],axis=0),color=tab10(rscount))
        axlat.scatter(np.mean(relstrat[stat],axis=1),yy[:,0],color=tab10(rscount))
        axlat.text(0.9,0.8-0.1*rscount,relstrat['label'],transform=axlat.transAxes,color=tab10(rscount),ha='right')
        
        if stat!='cv':
            print('stat is not cv')
            axlat.set_xscale('log')
            axlon.set_yscale('log')
        
        rscount+=1
    
    
    axlon.set_xlabel('Longitude')
    axlon.set_ylabel('particle wt')
    axlon.text(0.1,0.9,'{}) Long-avg {} particle wt'.format(atoz[rscount],stat),transform=axlon.transAxes)
    axlat.text(0.1,0.9,'{}) Lat-avg {} particle wt'.format(atoz[rscount+1],stat),transform=axlat.transAxes)
    axlat.set_xlabel('particle wt')
    axlat.set_ylabel('Latitude')
    # axlat.legend()
    
    plt.subplots_adjust(left=0.1,right=0.975,top=0.98,hspace=0.3,wspace=0.3)
    outfn = home + 'etools/plots/TV_vs_const_{}.png'.format(stat)
    plt.savefig(outfn)
    print(f'Saved to {outfn}')
    plt.close()
