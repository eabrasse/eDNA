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
from matplotlib import ticker
import matplotlib.dates as mdates
import pickle
import string
import efun
import pytz
from scipy.stats import linregress
import cmocean as cmo
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

tz = pytz.utc
pdt = pytz.timezone('America/Vancouver')

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('#526c30ff')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

blue_circle_color = '#37abc8ff'



plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'

atoz = string.ascii_lowercase

tab10 = plt.get_cmap('tab10',10)
Paired = plt.get_cmap('Paired',12)

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

# wrangle Feb samples, sigh.
transect_quant_fn = home+'data/All_Mod1_Molecular_Data - all_data.csv'
dfq = pd.read_csv(transect_quant_fn,sep=',',engine='python')
gq = dfq[dfq['campaign']=='Feb_2023']
ggq= gq[gq['task']=='transect']
# transect_metadata_fn = home+'data/MURI_Module1_qPCR_metadata - 02_samples.csv'
# dfm = pd.read_csv(transect_metadata_fn,sep=',',engine='python')
sample_dt = datetime(2023,2,3,13,0,tzinfo=pdt).astimezone(tz)


# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

#find June sampling particle tracks
data_fn = home+'data/Feb2023_hist_counts.p'

D = pickle.load(open(data_fn,'rb'))
tind = argnearest(D['dt_list0'],sample_dt)

moorrad_fn = home+'data/Feb2023_moor_radius_counts.p'
Dm = pickle.load(open(moorrad_fn,'rb'))
volume_list = [(np.pi*rad**2)*np.abs(Dm['zref']) for rad in Dm['radius_list']]
over_volume_list = [1/vol for vol in volume_list]
nrad = len(Dm['radius_list'])
nmoor = len(Dm['moor_list'])

cividis = plt.get_cmap('cividis',nrad)

D['bin_x_center'] = 0.5*(D['bin_x_edges'][:-1]+D['bin_x_edges'][1:])
D['bin_y_center'] = 0.5*(D['bin_y_edges'][:-1]+D['bin_y_edges'][1:])

fw,fh = efun.gen_plot_props()
fig = plt.figure(figsize=(fw*4,fh))
gs = GridSpec(1,4)
ax_samp =fig.add_subplot(gs[1])
ax_samp2 = ax_samp.twinx()
ax_age = fig.add_subplot(gs[2])
ax_deg_eff = fig.add_subplot(gs[3])
# ax_qq = plt.subplot(gs[3])

# # get grid from random file
grid_fn= home+ 'data/all_3d_hc_dolph_releases.p'

Dg = pickle.load(open(grid_fn,'rb'))

# gather some fields, for convenience
lonr = Dg['metadata']['lon'][:]
latr = Dg['metadata']['lat'][:]
maskr = Dg['metadata']['mask'][:]

xr,yr = efun.ll2xy(lonr,latr,lon0,lat0)

agescale = 1/3600

props = dict(boxstyle='round', facecolor='white',edgecolor='none',alpha=0.8)

k_degrade = -0.02

moor_count=0
particles = np.zeros((nmoor,nrad))
samples = np.zeros((nmoor))
distances = np.zeros((nmoor))
for moor in moor_list:
    
    gg = ggq[gq['reference']==moor['label']]
    sampledata = gg['dolphin_DNA_copy_uL'].values
    mean_sampledata = np.nanmean(sampledata)
    nsamples = len(sampledata)
    samples[moor_count] = mean_sampledata
    
    xm,ym = efun.ll2xy(moor['lon'],moor['lat'],lon0,lat0)
    moor['dist'] = np.sqrt(xm**2+ym**2)

    moor_rad = Dm['moor_list'][moor_count]
    
    conc_list = [moor_rad['count'][r]*over_volume_list[r] for r in range(nrad)]
    particles[moor_count,:] = conc_list
    
    distances[moor_count] = moor['dist']
    
    # ax_samp.plot(moor['dist'],samples[moor_count],marker='o',mfc=blue_circle_color,mec='k',linestyle='None',markersize=14,alpha=1.0,zorder=15)
    ax_samp.plot(moor['dist']*np.ones((nsamples)),sampledata,marker='o',mfc=blue_circle_color,mec='None',linestyle='None',markersize=10,alpha=0.5,zorder=5)
    # ax_samp.plot(moor['dist'],particles[moor_count],marker='*',mfc='yellow',mec='k',linestyle='None',markersize=10,alpha=1.0,zorder=15)
    
    for r in range(nrad):
        ax_age.text(0.3,.98-0.1*r,Dm['radius_list'][r],color=cividis(r),transform=ax_age.transAxes,va='top',ha='right')
    
        ax_age.plot(moor['dist'],agescale*Dm['particle_age_bins'][moor_count,r],marker='o',mfc=cividis(r),mec='None',linestyle='none',markersize=10,alpha=1.0)
        
        degradation_scale = np.exp(agescale*Dm['particle_age_bins'][moor_count,r]*k_degrade)
        
        ax_samp2.plot(moor['dist'],degradation_scale*conc_list[r],marker='*',mfc=cividis(r),mec='None',linestyle='None',markersize=10,alpha=1)
        
        # ax_deg_eff.plot(moor['dist'],conc_list[r]-degradation_scale*conc_list[r],marker='o',mfc=cividis(r),mec='None',linestyle='None',markersize=10,alpha=1)
        ax_deg_eff.plot(moor['dist'],degradation_scale,marker='o',mfc=cividis(r),mec='None',linestyle='None',markersize=10,alpha=1)
        
    # ax_qq.plot(particles[moor_count]*np.ones((nsamples)),sampledata,marker='o',mfc=blue_circle_color,mec='None',linestyle='None',markersize=10,alpha=0.5,zorder=5)
    
    moor_count+=1

fig2 = plt.figure(figsize=(fw*2,fh))
gs = GridSpec(1,2)
ax_r = fig2.add_subplot(gs[0])
ax_rmse = fig2.add_subplot(gs[1])

dinds = np.argsort(distances)
y0 = samples[dinds]
for r in range(nrad):
    
    x0 = particles[dinds,r]
    slope, intercept, r_value, p_value, std_err = linregress(x0,y0)
    ax_r.plot(Dm['radius_list'][r],r_value,marker='o',mfc=cividis(r),mec='None',linestyle='none',markersize=10,alpha=1.0)
    ax_rmse.plot(Dm['radius_list'][r],std_err,marker='o',mfc=cividis(r),mec='None',linestyle='none',markersize=10,alpha=1.0)
    
ax_r.set_xlabel('radius (m)')
ax_r.set_ylabel('R')
ax_r.text(0.1,0.9,'a) R value',transform=ax_r.transAxes)
ax_rmse.set_xlabel('radius (m)')
ax_rmse.set_ylabel('std err')
ax_rmse.text(0.1,0.9,'b) Standard Error',transform=ax_rmse.transAxes)

ax_samp.set_xlabel('Dist from pen (m)')
ax_samp.set_ylabel(r'DNA copies $\mathrm{\mu L}^{-1}$')

ax_samp.set_yscale('log')
ax_samp2.set_yscale('log')
ax_samp2.set_ylabel('Particle concentration')
ax_samp2.text(0.1,0.9,'INCLUDING DEGRADATION',transform=ax_samp2.transAxes)

ax_age.set_xlabel('Dist fom pen (m)')
ax_age.set_ylabel('Mean age of particles (hr)')

ax_deg_eff.set_xlabel('dist_from_pen')
ax_deg_eff.set_ylabel('degradation effect')
# MAP


ax = fig.add_subplot(gs[0])

ax.pcolormesh(xr,yr,maskr,cmap=cmap_mask,shading='nearest',zorder=100)
binx,biny = np.meshgrid(D['bin_x_center'],D['bin_y_center'])
p=ax.pcolormesh(binx,biny,D['hist_particles_bin'][tind,:],cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(),zorder=1,shading='nearest')

ax.set_xlabel('Dist from pen (m)')
ax.set_ylabel('Dist from pen (m)')
cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower right',bbox_transform=ax.transAxes,bbox_to_anchor=(-0.05,0,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
cb.set_label('Particles')

moor_count=0
for moor in moor_list:
    
        xmoor,ymoor = efun.ll2xy(moor['lon'],moor['lat'],lon0,lat0)
        ax.plot(xmoor,ymoor,marker='o',mec='k',mfc = blue_circle_color,markersize=10,zorder=200)

        moor_count+=1


ax.axis([-400,250,-500,750])

ax.set_aspect(1)
   
ax.text(0.1,0.9,'a)',transform=ax.transAxes,ha='left',va='top')
ax_samp.text(0.1,0.9,'b)',transform=ax_samp.transAxes,ha='left',va='top')


fig.subplots_adjust(right=0.98,left=0.02,bottom=0.15,top = 0.98,wspace=0.5)
# fig2.subplots_adjust(right=0.98,left=0.2,bottom=0.2,top = 0.98)
# outfn = home+'plots/Feb2023_moorings_fit.png'
# plt.savefig(outfn)
# print(f'Saved to {outfn}')
plt.show(block=False)
plt.pause(0.1)
# plt.close()


