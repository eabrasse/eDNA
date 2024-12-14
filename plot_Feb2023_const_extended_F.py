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
ggq['copies_per_mLseawater'] = ggq.apply(lambda row: row.dolphin_DNA_copy_uL*row.volume_eluted/(row.volume_L*1000),axis=1)
sample_dt = datetime(2023,2,3,13,0,tzinfo=pdt).astimezone(tz)


# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

#find June sampling particle tracks
data_fn = home+'data/Feb2023_hist_counts_dev2.p'

D = pickle.load(open(data_fn,'rb'))
tind = argnearest(D['dt_list0'],sample_dt)

moorrad_fn = home+'data/Feb2023_moor_radius_counts_dev2.p'
Dm = pickle.load(open(moorrad_fn,'rb'))
volume_list = [(np.pi*rad**2)*np.abs(Dm['zref']) for rad in Dm['radius_list']]
over_volume_list = [1/vol for vol in volume_list]
nrad = len(Dm['radius_list'])
nmoor = len(Dm['moor_list'])

D['bin_x_center'] = 0.5*(D['bin_x_edges'][:-1]+D['bin_x_edges'][1:])
D['bin_y_center'] = 0.5*(D['bin_y_edges'][:-1]+D['bin_y_edges'][1:])

fw,fh = efun.gen_plot_props()
fig = plt.figure(figsize=(fw*2,fh*1.5))
gs = GridSpec(5,2)
ax_samp = plt.subplot(gs[:3,:])
ax_samp2 = ax_samp.twinx()
ax_age = plt.subplot(gs[-2:,0])
ax_qq = plt.subplot(gs[-2:,1])

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

moor_count=0
particles = np.zeros((nmoor))
samples = np.zeros((nmoor))
distances = np.zeros((nmoor))
age_meds = np.zeros((nmoor))
r=0
for moor in moor_list:
    
    gg = ggq[gq['reference']==moor['label']]
    sampledata = gg['copies_per_mLseawater'].values
    mean_sampledata = np.nanmean(sampledata)
    nsamples = len(sampledata)
    samples[moor_count] = mean_sampledata
    
    xm,ym = efun.ll2xy(moor['lon'],moor['lat'],lon0,lat0)
    moor['dist'] = np.sqrt(xm**2+ym**2)
    
    # i = np.argmin(np.abs(D['bin_x_edges']-xm))
    # j = np.argmin(np.abs(D['bin_y_edges']-ym))
    # particles[moor_count] = D['hist_particles_bin'][tind,j,i]
    
    moor_rad = Dm['moor_list'][moor_count]
    
    conc_list = moor_rad['count'][r]*over_volume_list[r]
    particles[moor_count] = conc_list

    distances[moor_count] = moor['dist']
    
    ax_samp.plot(moor['dist'],samples[moor_count],marker='o',mfc=blue_circle_color,mec='k',linestyle='None',markersize=14,alpha=1.0,zorder=15)
    ax_samp.plot(moor['dist']*np.ones((nsamples)),sampledata,marker='o',mfc=blue_circle_color,mec='None',linestyle='None',markersize=10,alpha=0.5,zorder=10)
    ax_samp2.plot(moor['dist'],conc_list,marker='*',mfc='orange',mec='k',linestyle='None',markersize=14,alpha=1.0,zorder=15)
    
    age_list = Dm['particle_age_lists'][r][moor_count]
    age_list_scaled = [agescale*age for age in age_list]
    
    age_meds[moor_count] = np.median(age_list_scaled)
    
    ax_age.scatter([moor['dist']]*len(age_list_scaled),age_list_scaled,marker='*',c='orange',s=40,alpha=0.5)
    ax_age.boxplot(age_list_scaled,whis=(0,100),positions=[moor['dist']],vert=True, widths=100,medianprops=dict(linestyle='-', linewidth=2, color='k'))
    
    
    moor_count+=1

dind = np.argsort(distances)
x0 = distances[dind]
for ax,var,col in [ax_samp,samples,blue_circle_color],[ax_samp2,particles,'orange']:
    y0 = np.log(var[dind])
    slope, intercept, r_value, p_value, std_err = linregress(x0,y0)
    xx = np.linspace(x0.min(),x0.max(),100)
    yy = np.exp(xx*slope)*np.exp(intercept)
    ax.plot(xx,yy,marker='none',linestyle='dashed',color=col,zorder=2)
    print(f'slope={slope}')
# y0 = samples[pinds]
# # p = np.polyfit(x0,y0,1)
# slope, intercept, r_value, p_value, std_err = linregress(x0,y0)
# xx = np.linspace(x0.min(),x0.max(),100)
# yy = xx*slope+intercept
# ax_samp.plot(distances[dind],samples[dind],marker='none',linestyle='solid',color=blue_circle_color,zorder=10)
# ax_samp2.plot(distances[dind],particles[dind],marker='none',linestyle='solid',color='orange',zorder=10)
aind = np.argsort(age_meds)
ax_qq.plot(age_meds[aind],samples[aind]/samples.max(),marker='o',mfc=blue_circle_color,mec='k',linestyle='dashed',color=blue_circle_color,markersize=14,alpha=1.0,zorder=15)
ax_qq.plot(age_meds[aind],particles[aind]/particles.max(),marker='*',mfc='orange',mec='k',linestyle='dashed',color='orange',markersize=14,alpha=1.0,zorder=15)
xx = np.linspace(age_meds.min(),age_meds.max(),100)
ydegrade = np.exp(-0.02*xx)
ax_qq.plot(xx,ydegrade,color='gray',linestyle='dashed',marker='None')


posdata_fn = home+'data/Feb2023_surfaceparticlepositionstats.p'
Dpos = pickle.load(open(posdata_fn,'rb'))
t_list = np.arange(len(Dpos.keys()))
all_release_densities = {}
for t in t_list:
    all_release_densities[t] = []
for key in Dpos.keys():
    rel = Dpos[key]
    
    density = 1/(2*rel['avg_dist_from_COM']**2)
    ax_qq.plot(agescale*rel['t'],density/density[0],color=tab10(2),alpha=0.1,zorder=1)
    
    for t in agescale*rel['t']:
        all_release_densities[t].append(density[int(t)])

# if I expect them to decay exponentially...
mean_release_density = [np.mean(all_release_densities[t]) for t in t_list]
frac_release_density = mean_release_density/mean_release_density[0]
ax_qq.plot(t_list,frac_release_density,color=tab10(2),alpha=1,zorder=5,linestyle='dashed',marker='None')

# pinds = np.argsort(particles)
# x0 = particles[pinds]
# y0 = samples[pinds]
# # p = np.polyfit(x0,y0,1)
# slope, intercept, r_value, p_value, std_err = linregress(x0,y0)
# xx = np.linspace(x0.min(),x0.max(),100)
# yy = xx*slope+intercept
# moor_count=0
# for moor in moor_list:
#     ax_samp.plot(moor['dist'],slope*particles[moor_count],marker='*',mfc='orange',mec='k',linestyle='None',markersize=14,alpha=1.0,zorder=15)
#     moor_count+=1
ax_samp.set_xlabel('Dist (m)')
ax_samp.set_ylabel(r'DNA copies $\mathrm{mL}^{-1}$',color='#37abc8ff')
ax_samp.tick_params(axis='y',which='both',labelcolor='#37abc8ff')

ax_samp.set_yscale('log')
ax_samp2.set_yscale('log')
ax_samp2.set_ylabel(r'Particles $\mathrm{m}^{-3}$',color = 'orange')
ax_samp2.tick_params(axis='y',which='both',labelcolor='orange')

ax_age.set_xlabel('Dist (m)')
ax_age.set_ylabel('Age (hr)')

ax_age.set_xlim([-50,650])
xtick = range(0,601,100)
ax_age.set_xticks(xtick)
ax_age.set_xticklabels([f'{xt:}' for xt in xtick])
# ax_age.set_xticklabels(xtl)
ax_qq.set_xlabel('Age (hr)')
# ax_qq.set_ylabel('Fraction of concentration remaining')
# ax_qq.text(0.2,0.8,'Fraction remaining',transform=ax_qq.transAxes)
ax_qq.set_yscale('log')
ax_qq.set_ylabel('C/C0')
ax_qq.set_xlim([-1,age_meds.max()+3])

# MAP


ax = inset_axes(ax_samp,height='100%',width='100%',loc='upper right',\
    bbox_to_anchor=(0.6, .5, .6, .5),bbox_transform=ax_samp.transAxes)

ax.pcolormesh(xr,yr,maskr,cmap=cmap_mask,shading='nearest',zorder=100)
binx,biny = np.meshgrid(D['bin_x_center'],D['bin_y_center'])
p=ax.pcolormesh(binx,biny,D['hist_particles_bin'][tind,:],cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(),zorder=1,shading='nearest')
# pl=ax.contour(xbins,ybins,slope*moor_dict['const_DNA_bin'][sample_dt_ind,:].T,colors=['k'],linewidths=[0.5],locator=ticker.LogLocator(),zorder=2)
# ax.contour(xgrid,ygrid,maskr,levels=[0.5],colors=['k'],linewidths=[1.5],zorder=150)
ax.set_xlabel('Dist (m)',fontsize=8)
ax.set_ylabel('Dist (m)',fontsize=8)
ax.tick_params(axis='both',which='both',labelsize=8,size=2)
cbaxes = inset_axes(ax, width="50%", height="4%", loc='upper left',bbox_transform=ax.transAxes,bbox_to_anchor=(0.05,-0.05,1,1))
cbaxes.tick_params(axis='both',which='both',labelsize=8,size=2)
cb = fig.colorbar(p, cax=cbaxes,orientation='horizontal')
# cb.set_label(r'particles $\mathrm{m}^{-3}$')
# cb.set_label('Particles',fontsize=8)
# cb.add_lines(pl)
moor_count=0
for moor in moor_list:
    
        xmoor,ymoor = efun.ll2xy(moor['lon'],moor['lat'],lon0,lat0)
        ax.plot(xmoor,ymoor,marker='o',mec='k',mfc = blue_circle_color,markersize=5,zorder=200)

        moor_count+=1

ax.axis([-400,250,-500,750])

ax.set_aspect(1)

# ax_samp.text(0.1,0.9,'b)',transform=ax_samp.transAxes,ha='left',va='top')
yl0 = ax_samp.get_ylim()
yl1 = ax_samp2.get_ylim()
ymin = np.min([yl0[0],yl1[0]])
ymax = np.max([yl0[1],yl1[1]])
ylim = [ymin,ymax]
for ax in ax_samp, ax_samp2:
    ax.set_ylim(ylim)

fig.subplots_adjust(right=0.85,left=0.15,bottom=0.1,top = 0.98,wspace=0.4,hspace=0.9)
# fig2.subplots_adjust(right=0.98,left=0.2,bottom=0.2,top = 0.98)
outfn = home+'plots/Feb2023_moorings_fit_radius_boxplot_dev.png'
plt.savefig(outfn,format='png',dpi=600)
# print(f'Saved to {outfn}')
# plt.show(block=False)
# plt.pause(0.1)
plt.close()


