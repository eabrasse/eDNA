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

# load mooring results
moor_fn = home+'LO_data/eDNA/Feb2023_moorings.p'

moor_dict = pickle.load(open(moor_fn,'rb'))
const= {}
const['particle_bin'] = moor_dict['const_particle_bin']
const['label'] = 'Constant'
TV = {}
TV['particle_bin'] = moor_dict['TV_particle_bin']
TV['label'] = 'Time-varying'
relstrat_list = [const,TV]
for relstrat in relstrat_list:
    relstrat['mean'] = np.mean(relstrat['particle_bin'],axis=0)
    relstrat['std'] = np.std(relstrat['particle_bin'],axis=0)
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
colmap = cmo.cm.matter
normal = matplotlib.colors.LogNorm()
for stat in statfig_list:
    
    plt.close('all')
    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(2,len(relstrat_list)+1)
    axlon = fig.add_subplot(gs[1,len(relstrat_list)])
    axlat = fig.add_subplot(gs[0,len(relstrat_list)])
    
    rscount=0
    for relstrat in relstrat_list:
        axmap = fig.add_subplot(gs[:,rscount])
        axmap.contour(lonp,latp,maskr,levels=[0.5],colors='k',linewidths=1,linestyles='solid')
        
        p=axmap.pcolormesh(xx,yy,relstrat[stat],cmap = colmap,norm=normal) 
        cbaxes = inset_axes(axmap, width="4%", height="60%", loc='center right',bbox_transform=axemap.transAxes,bbox_to_anchor=(0.15,0.,1,1))
        cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
        cb.set_label('{} particle wt'.format(stat))
        
        axmap.text(0.1,0.9,'{}) {} {}'.format(atoz[rscount],relstrat['label'],stat),transform=axmap.transAxes)
        axmap.set_xlabel('Longitude')
        axmap.set_ylabel('Latitude')
        
        axlon.plot(xx,relstrat[stat][0,:],color=tab10(rscount),label=relstrat['label'])
        axlat.plot(relstrat[stat][:,0],yy,color=tab10(rscount),label=relstrat['label'])
        
        rscount+=1
        
    axlon.set_xlabel('Longitude')
    axlon.set_ylabel('Long-avg {} particle wt'.format(stat))
    axlat.set_xlabel('Lat-avg {} particle wt'.format(stat))
    axlat.set_ylabel('Latitude')
    axlat.legend()
    
    outfn = home + 'etools/plots/TV_vs_const_{}.png'.format(stat)
    plt.savefig(outfn)
    print(f'Saved to {outfn}')
    plt.close()
