#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot map of COAWST/CSIDE/SD Bight model with landmarks for context
"""

# setup
import os
import sys
import netCDF4 as nc
import efun
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import string
from PIL import Image
s10 = plt.get_cmap('Set1')

plt.close('all')
fs_small=10
fs_big = 12
atoz=string.ascii_lowercase

def add_features(ax,axes):
    # use cartopy features to plot landmarks on map
    ax.set_extent(axes, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    # ax.add_feature(cfeature.RIVERS)
    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

    SOURCE = 'Natural Earth'
    LICENSE = 'public domain'

    ax.add_feature(states_provinces, edgecolor='gray')

def plot_domain(lonr,latr,ax,col):
    # outline different model level domains
    # note: this program was more useful when I was plotting all 4 domains.
    # eventually, I only decided to plot the relevant LV4 domain
    lw = 1
    ls = 'solid'
    ax.plot(lonr[0,:],latr[0,:],color=col,lw=lw,linestyle=ls)
    ax.plot(lonr[-1,:],latr[-1,:],color=col,lw=lw,linestyle=ls)
    ax.plot(lonr[:,0],latr[:,0],color=col,lw=lw,linestyle=ls)
    ax.plot(lonr[:,-1],latr[:,-1],color=col,lw=lw,linestyle=ls)

props = dict(boxstyle='round', fc='white',ec='None',alpha=1.0)

# load in model domains
home = '/Users/elizabethbrasseale/Projects/eDNA/'
grid_dir = home+'data/'
gridfn = grid_dir + 'grid.nc'
dsg = nc.Dataset(gridfn)

# load background map
figname = home+'data/bangor_gmaps.png'
img_extent = (-122.754167, -122.695833, 47.718611, 47.788056)

# define axis ranges
# note: I previously included an inset axis showing all 4 domains. 
axes_big = [-129,-122,45,53]
# axes_med = [-123.5,-122,46.5,49]
axes_med = [-123,-122.25,47,48.5]
# axes_small = [-117.23,-117.09,32.43,32.7]

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

# generate figure and name axis handles
fw,fh = efun.gen_plot_props(fs_big=fs_big,fs_small=fs_small)
fig = plt.figure(figsize=(1.75*fw,fh))
ax0 = fig.add_subplot(1,2,1,projection=ccrs.PlateCarree())
ax1 = fig.add_subplot(1,2,2,projection=ccrs.PlateCarree())
# ax2 = fig.add_subplot(1,3,3,projection=ccrs.PlateCarree())
# axs = [ax0,ax1]

# ax0 = plt.axes(projection=ccrs.PlateCarree())
count=0
for ax,axes in [[ax0,axes_big],[ax1,axes_med]]:#,[ax2,img_extent]]:
    add_features(ax=ax,axes=axes)
    ax.plot(lon0,lat0,mfc='magenta',mec='k',marker='*',markeredgewidth=0.75,markersize=15)
    ax.text(0.05,0.05,'{})'.format(atoz[count]),ha='left',va='bottom',transform=ax.transAxes)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.tick_params(length=3)
    count+=1
# add_features(ax=ax0,axes=axes_big)
# ax0.plot(lon0,lat0,mfc=s10(4),mec='k',marker='*',markersize=18)
# ax0.text(0.05,0.05,'a)',ha='left',va='bottom',transform=ax0.transAxes)
# ax0.set_xlabel('Longitude',fontsize=10)
# ax0.set_ylabel('Latitude',fontsize=10)
ax0.set_xticks([-127,-125,-123])
ax0.set_yticks([45,47,49,51,53])
# ax0.tick_params(length=3)
ax1.set_xticks([-123,-122.5])
ax1.set_yticks([47,47.5,48,48.5])

# with Image.open(figname) as img:
#     ax2.imshow(img, extent=img_extent,transform=ccrs.PlateCarree())
#     ax2.plot(lon0,lat0,mfc='magenta',mec='k',markeredgewidth=0.5,marker='*',markersize=16)
#     ax2.set_xlabel('Longitude')
#     ax2.set_ylabel('Latitude')
# ax2.set_xticks([-122.75,-122.725,-122.7])
# ax2.set_yticks([47.73,47.75,47.77])
# ax2.tick_params(length=3)
# ax2.text(0.05,0.05,'{})'.format(atoz[count]),ha='left',va='bottom',transform=ax2.transAxes,color='w')
# gl=ax2.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False,ylocs=[47.73,47.75,47.77],xlocs=[-122.75,-122.725,-122.7])
# gl.top_labels = False
# gl.right_labels = False

# show or save plot
plt.subplots_adjust(bottom=0.13,left=0.1,right=0.98,wspace=0.1)
# plt.show(block=False)
# plt.pause(0.1)
outfn = home+'plots/contextmap2.png'
plt.savefig(outfn,format='png',transparent=True)