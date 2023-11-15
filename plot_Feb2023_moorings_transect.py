import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from PIL import Image
import numpy as np
from datetime import datetime, timedelta
import pytz
import netCDF4 as nc
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.dates as mdates
import pickle
import string
import efun
import pytz

tz = pytz.utc
pdt = pytz.timezone('America/Vancouver')

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'

atoz = string.ascii_lowercase

tab10 = plt.get_cmap('tab10',10)
Paired = plt.get_cmap('Paired',12)


# wrangle Feb samples, sigh.
transect_quant_fn = home+'data/MURI_Module1_datasheets_07_transect_Feb2023.csv'
dfq = pd.read_csv(transect_quant_fn,sep=',',engine='python')
gq = dfq[dfq['sample_type']=='eDNA']
mll2ref = {'Marginal Pier':'Marginal Pier ',
    'near Marginal Pier':'Half way between Marginal Pier and bridge',
    'near Delta Pier':'Near bridge',
    'Delta Pier':'Outside dolphin heated pens from boat'}
# transect_metadata_fn = home+'data/MURI_Module1_qPCR_metadata - 02_samples.csv'
# dfm = pd.read_csv(transect_metadata_fn,sep=',',engine='python')
sample_dt = datetime(2023,2,3,13,0,tzinfo=pdt).astimezone(tz)


# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

#find June sampling particle tracks
moor_fn = home+'data/Feb2023_DNA_moorings.p'

moor_dict = pickle.load(open(moor_fn,'rb'))
moor_list = moor_dict['moor_list']
dt_list = moor_dict['dt_list']
nmoor = len(moor_list)
rainbow = plt.get_cmap('rainbow',nmoor-1)
moor_lat_list = [moor['lat'] for moor in moor_list]
ns_inds = np.argsort(moor_lat_list)[::-1]

#get mean concentration
data_fn = home+'data/03_ESP1_Feb2023_hourly.csv'
df = pd.read_csv(data_fn,sep=',',engine='python')
DNA_mean = np.mean(df.PB_quantity_mean)


fig = plt.figure(figsize=(12, 8))
gs = GridSpec(nmoor-1+1,3)

ax = plt.subplot(gs[:,0])

# # get grid from random file
grid_fn= home+ 'data/all_3d_hc_dolph_releases.p'

D = pickle.load(open(grid_fn,'rb'))

# gather some fields, for convenience
lonr = D['metadata']['lon'][:]
latr = D['metadata']['lat'][:]
maskr = D['metadata']['mask'][:]


props = dict(boxstyle='round', facecolor='white',edgecolor='none',alpha=0.8)
#calculate release volume
# key_list = D.keys()
# release_list = [key for key in key_list if key[:3]=='rel']
# example_release = D[release_list[0]]
# lon = example_release['lon'][:]
# lat = example_release['lat'][:]
# z = example_release['z'][:]
# get number of particles
# particle_rel = lon.shape[0]*41
# # derive the area over which they were released
# xp0,yp0 = efun.ll2xy(lon,lat,lon0,lat0)
# dxrel = np.abs(xp0.max()-xp0.min())
# dyrel = np.abs(yp0.max()-yp0.min())
# dzrel = 1 #all released at surface, soooo....
# vol_rel = dxrel*dyrel*dzrel
#
# particle_conc_rel = particle_rel/vol_rel

#calculate bin volumes
# set domain limits
# pad = .01
# # plot region around delta pier, using release point as indicator
# aa = [lon0-pad, lon0+pad,
#    lat0-pad, lat0+pad]
# aax,aay =  efun.ll2xy(np.array(aa[:2]),np.array(aa[2:]),lon0,lat0)
# aaxy = [aax[0],aax[1],aay[0],aay[1]]
# zref = -2
# nbins = 100
# bin_x_edges=np.linspace(aax[0], aax[1],nbins+1)
# bin_y_edges=np.linspace(aay[0], aay[1],nbins+1)
# xx, yy = np.meshgrid(bin_x_edges[:-1]+0.5*(bin_x_edges[1]-bin_x_edges[0]),bin_y_edges[:-1]+0.5*(bin_y_edges[1]-bin_y_edges[0]))
# dxbin = bin_x_edges[1]-bin_x_edges[0]
# dybin = bin_y_edges[1]-bin_y_edges[0]
# dzbin = np.abs(zref)
# vol_bin = dxbin*dybin*dzbin

# MAP
xgrid,ygrid = efun.ll2xy(lonr,latr,lon0,lat0)
ax.pcolormesh(xgrid,ygrid,maskr,cmap=cmap_mask,shading='nearest',zorder=0)
ax.contour(xgrid,ygrid,maskr,levels=[0.5],colors=['k'],linewidths=[1.5])
ax.set_xlabel('Dist from pen (m)')
ax.set_ylabel('Dist from pen (m)')

moor_count=0
mooring_axes = []

for ind in ns_inds:

    moor = moor_list[ind]
    if moor['label'] in mll2ref.keys():
        #calculate DNA at mooring
        # const_DNA_particle_conc_bin = DNA_mean*moor['const_particle_bin']/vol_bin
        # const_DNA_bin = const_DNA_particle_conc_bin/particle_conc_rel
        # TV_DNA_particle_conc_bin = moor['TV_particle_bin']/vol_bin #note: already scaled by DNA
        # TV_DNA_bin = TV_DNA_particle_conc_bin/particle_conc_rel
    
    
        x,y = efun.ll2xy(moor['lon'],moor['lat'],lon0,lat0)
        ax.plot(x,y,marker='d',mec='k',mfc = rainbow(moor_count),markersize=10)

        axm = plt.subplot(gs[moor_count,1:])

        # moor['const_particle_bin'][moor['const_particle_bin']==0] = np.nan
        axm.plot(dt_list,moor['const_DNA_bin'],linestyle='dashed',color=rainbow(moor_count))
        # moor['TV_particle_bin'][moor['TV_particle_bin']==0] = np.nan
        axm.plot(dt_list,moor['TV_DNA_bin'],linestyle='solid',color=rainbow(moor_count))

    
        axm.text(0.05,0.8,'{}) {}'.format(atoz[moor_count],moor['label']),color='k',transform=axm.transAxes,ha='left',va='top')
        axm.grid()

        # get universal ylims
        if moor_count==0:
            ymax = max([np.nanmax(moor['const_DNA_bin']),np.nanmax(moor['TV_DNA_bin'])])
            ymin = min([np.nanmin(moor['const_DNA_bin'][moor['const_DNA_bin']>0]),np.nanmin(moor['TV_DNA_bin'][moor['TV_DNA_bin']>0])])
        ymax = max([np.nanmax(moor['const_DNA_bin']),np.nanmax(moor['TV_DNA_bin']),ymax])
        ymin = min([np.nanmin(moor['const_DNA_bin'][moor['const_DNA_bin']>0]),np.nanmin(moor['TV_DNA_bin'][moor['TV_DNA_bin']>0]),ymin])

        # format date/time axis
        # if moor_count<(nmoor-2):
        axm.set_xlabel('')
        axm.set_xticklabels(['' for xtl in axm.get_xticklabels()])
        # else:
            # axm.set_xlabel('Date & time')
            # axm.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d %H:%m"))
            # plt.setp( axm.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')

        mooring_axes.append(axm)
    
    
        reference=mll2ref[moor['label']]
        gg = gq[gq['reference']==reference]
        sampledata = gg['PB_quantity'].values
        ndatapoints = len(sampledata)
        datatime = [sample_dt]*ndatapoints
        axm.plot(datatime,sampledata,marker='o',mfc=rainbow(moor_count),mec='k',linestyle='None',markersize=10,alpha=0.5)
        ymax = max([np.nanmax(sampledata),ymax])
        ymin = min([np.nanmin(sampledata),ymin])
        
        mean_sampledata = np.nanmean(sampledata)
        axm.axhline(y=mean_sampledata,linestyle='dotted',color=rainbow(moor_count))

    # if moor['label']=='Dolphin pen':
    #     axm.axhline(y=DNA_mean,linestyle='dotted',color=rainbow(moor_count))
        if moor_count==0:
            axm.text(0.8,0.95,'Solid = time-varying\nDashed = constant\nDots = samples\nDotted line = sample mean',transform=axm.transAxes,ha='right',va='top',bbox=props,fontsize=10)

        moor_count+=1

axp = plt.subplot(gs[-1,1:])
axp.fill_between(dt_list,np.zeros((len(dt_list))),moor_dict['active_particles'],color='gray')
axp.set_xlabel('Date & time (UTC)')
axp.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d %H:%m"))
plt.setp( axp.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
axp.text(0.05,0.8,'{}) {}'.format(atoz[moor_count],'Active particles'),color='k',transform=axp.transAxes,ha='left',va='top')
axp.grid()

# ax.axis([-122.735,-122.725,47.74,47.75])
lonaxes = np.array([-122.735,-122.725]);lataxes=np.array([47.74,47.75])
xaxes,yaxes = efun.ll2xy(lonaxes,lataxes,lon0,lat0)
ax.axis([xaxes[0],xaxes[1],yaxes[0],yaxes[1]])
ax.set_aspect(1)
   

for axm in mooring_axes:
    axm.set_ylim([ymin,ymax])
    axm.set_yscale('log')
    axm.set_ylabel('DNA conc\n'+r'(copies $\mathrm{\mu L}^{-1}}$)')
    # print('ymin = {}, ymax = {}'.format(ymin, ymax))


fig.subplots_adjust(right=0.98,left=0.1,bottom=0.1,top = 0.98,wspace=0.5)
# outfn = home+'plots/Feb2023_moorings_obs.png'
# plt.savefig(outfn)
# print(f'Saved to {outfn}')
plt.show(block=False)
plt.pause(0.1)
# plt.close()


