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

D['bin_x_center'] = 0.5*(D['bin_x_edges'][:-1]+D['bin_x_edges'][1:])
D['bin_y_center'] = 0.5*(D['bin_y_edges'][:-1]+D['bin_y_edges'][1:])

fw,fh = efun.gen_plot_props()
fig = plt.figure(figsize=(fw*4,fh))
gs = GridSpec(1,4)
ax_samp = plt.subplot(gs[1])
ax_age = plt.subplot(gs[2])
ax_qq = plt.subplot(gs[3])

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
for moor in moor_list:
    
    gg = ggq[gq['reference']==moor['label']]
    sampledata = gg['dolphin_DNA_copy_uL'].values
    mean_sampledata = np.nanmean(sampledata)
    nsamples = len(sampledata)
    samples[moor_count] = mean_sampledata
    
    xm,ym = efun.ll2xy(moor['lon'],moor['lat'],lon0,lat0)
    moor['dist'] = np.sqrt(xm**2+ym**2)
    
    i = np.argmin(np.abs(D['bin_x_edges']-xm))
    j = np.argmin(np.abs(D['bin_y_edges']-ym))
    particles[moor_count] = D['hist_particles_bin'][tind,j,i]

    
    
    
    
    distances[moor_count] = moor['dist']
    
    ax_samp.plot(moor['dist'],samples[moor_count],marker='o',mfc=blue_circle_color,mec='k',linestyle='None',markersize=14,alpha=1.0,zorder=15)
    ax_samp.plot(moor['dist']*np.ones((nsamples)),sampledata,marker='o',mfc=blue_circle_color,mec='None',linestyle='None',markersize=10,alpha=0.5,zorder=5)
    # ax_samp.plot(moor['dist'],particles[moor_count],marker='*',mfc='yellow',mec='k',linestyle='None',markersize=10,alpha=1.0,zorder=15)
    
    ax_age.plot(moor['dist'],agescale*D['particle_age_bins'][tind,j,i],marker='o',mfc=tab10(1),mec='None',markersize=10)
    ax_qq.plot(particles[moor_count]*np.ones((nsamples)),sampledata,marker='o',mfc=blue_circle_color,mec='None',linestyle='None',markersize=10,alpha=0.5,zorder=5)
    
    moor_count+=1
        
pinds = np.argsort(particles)
x0 = particles[pinds]
y0 = samples[pinds]
# p = np.polyfit(x0,y0,1)
slope, intercept, r_value, p_value, std_err = linregress(x0,y0)
xx = np.linspace(x0.min(),x0.max(),100)
yy = xx*slope+intercept
moor_count=0
for moor in moor_list:
    ax_samp.plot(moor['dist'],slope*particles[moor_count],marker='*',mfc='yellow',mec='k',linestyle='None',markersize=14,alpha=1.0,zorder=15)
    moor_count+=1
# ax_samp.plot(xx,yy,linestyle='dashed',color='k',zorder=10)
# ax_samp.text(0.5,0.5,'best fit line = {:0.0f}x+{:0.1f}'.format(slope,intercept),transform=ax_samp.transAxes,ha='center',va='center')
ax_samp.set_xlabel('Dist from pen (m)')
ax_samp.set_ylabel(r'DNA copies $\mathrm{\mu L}^{-1}$')
# ax_samp.set_xscale('log')
ax_samp.set_yscale('log')
# ax_samp2.set_yscale('log')
# ax_samp.text(0.4,0.65,f'Linear fit\nslope = {slope:2.0f}',transform=ax_samp.transAxes,rotation=30)
# ax_samp.text(0.4,0.4,r'$\mathrm{R}^{2}$ '+'= {:.4}'.format(r_value**2)+'\n'+'p={:.2E}'.format(p_value),transform=ax_samp.transAxes)
# ax_samp2.set_ylabel(r'particles $\mathrm{m}^{-3}$')

ax_age.set_xlabel('Dist fom pen (m)')
ax_age.set_ylabel('Mean age of particles (hr)')
ax_qq.set_xlabel('Particles')
ax_qq.set_ylabel(r'Samples (copies $\mathrm{\mu L}^{-1}$)')
ax_qq.set_yscale('log')
ax_qq.set_xscale('log')

# MAP


ax = plt.subplot(gs[0])
# xgrid,ygrid = efun.ll2xy(lonr,latr,lon0,lat0)
# #~~~~~~~~~~~~
# #reprint code to calculate bins because I forgot to output them
# # set domain limits
# pad = .01
# # plot region around delta pier, using release point as indicator
# aa = [lon0-pad, lon0+pad,
#    lat0-pad, lat0+pad]
# aax,aay =  efun.ll2xy(np.array(aa[:2]),np.array(aa[2:]),lon0,lat0)
# aaxy = [aax[0],aax[1],aay[0],aay[1]]
# nbins = 100
# bin_x_edges=np.linspace(aax[0], aax[1],nbins+1)
# xbins = 0.5*(bin_x_edges[:-1]+bin_x_edges[1:])
# bin_y_edges=np.linspace(aay[0], aay[1],nbins+1)
# ybins = 0.5*(bin_y_edges[:-1]+bin_y_edges[1:])
# #~~~~~~~~~~~~~
ax.pcolormesh(xr,yr,maskr,cmap=cmap_mask,shading='nearest',zorder=100)
binx,biny = np.meshgrid(D['bin_x_center'],D['bin_y_center'])
p=ax.pcolormesh(binx,biny,D['hist_particles_bin'][tind,:],cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(),zorder=1,shading='nearest')
# pl=ax.contour(xbins,ybins,slope*moor_dict['const_DNA_bin'][sample_dt_ind,:].T,colors=['k'],linewidths=[0.5],locator=ticker.LogLocator(),zorder=2)
# ax.contour(xgrid,ygrid,maskr,levels=[0.5],colors=['k'],linewidths=[1.5],zorder=150)
ax.set_xlabel('Dist from pen (m)')
ax.set_ylabel('Dist from pen (m)')
cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower right',bbox_transform=ax.transAxes,bbox_to_anchor=(-0.05,0,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
# cb.set_label(r'particles $\mathrm{m}^{-3}$')
cb.set_label('Particles')
# cb.add_lines(pl)
moor_count=0
for moor in moor_list:
    
        xmoor,ymoor = efun.ll2xy(moor['lon'],moor['lat'],lon0,lat0)
        ax.plot(xmoor,ymoor,marker='o',mec='k',mfc = blue_circle_color,markersize=10,zorder=200)

        # # get universal ylims
        # if moor_count==0:
        #     ymax = np.nanmax(slope*moor['const_DNA_bin'])
        #     ymin = np.nanmin(slope*moor['const_DNA_bin'][moor['const_DNA_bin']>0])
        # ymax = max([np.nanmax(slope*moor['const_DNA_bin']),ymax])
        # ymin = min([np.nanmin(slope*moor['const_DNA_bin'][moor['const_DNA_bin']>0]),ymin])

    
        # reference=mll2ref[moor['label']]
        # gg = gq[gq['reference']==reference]
        # sampledata = gg['dolphin_DNA_copy_uL'].values
        # ndatapoints = len(sampledata)
        # datatime = [sample_dt]*ndatapoints
        # # axm.plot(datatime,sampledata,marker='o',mfc=rainbow(moor_count),mec='k',linestyle='None',markersize=10,alpha=0.5)
        # ymax = max([np.nanmax(sampledata),ymax])
        # ymin = min([np.nanmin(sampledata),ymin])
        #
        # mean_sampledata = np.nanmean(sampledata)
        # axm.axhline(y=mean_sampledata,linestyle='dotted',color=rainbow(moor_count))
        

    # if moor['label']=='Dolphin pen':
    #     axm.axhline(y=DNA_mean,linestyle='dotted',color=rainbow(moor_count))
        # if moor_count==0:
            # axm.text(0.8,0.95,'Solid = time-varying\nDashed = constant\nDots = samples\nDotted line = sample mean',transform=axm.transAxes,ha='right',va='top',bbox=props,fontsize=10)

        moor_count+=1


# axp = plt.subplot(gs[-1,1:])
# axp.fill_between(dt_list,np.zeros((len(dt_list))),moor_dict['active_particles'],color='gray')
# axp.set_xlabel('Date & time (UTC)')
# axp.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d %H:%m"))
# plt.setp( axp.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
# axp.text(0.05,0.8,'{}) {}'.format(atoz[moor_count],'Active particles'),color='k',transform=axp.transAxes,ha='left',va='top')
# axp.grid()

# ax.axis([D['bin_x_edges'].min(),D['bin_x_edges'].max(),D['bin_y_edges'].min(),D['bin_y_edges'].max()])
ax.axis([-400,250,-500,750])
# lonaxes = np.array([-122.735,-122.725]);lataxes=np.array([47.74,47.75])
# xaxes,yaxes = efun.ll2xy(lonaxes,lataxes,lon0,lat0)
# pad = 300
# ax.axis([xaxes[0]-pad,xaxes[1]+pad,yaxes[0]-pad,yaxes[1]+pad])
ax.set_aspect(1)
   
ax.text(0.1,0.9,'a)',transform=ax.transAxes,ha='left',va='top')
ax_samp.text(0.1,0.9,'b)',transform=ax_samp.transAxes,ha='left',va='top')
# for axm in mooring_axes:
#     axm.set_ylim([ymin,ymax])
#     axm.set_yscale('log')
#     axm.set_ylabel('DNA conc\n'+r'(copies $\mathrm{\mu L}^{-1}}$)')
    # print('ymin = {}, ymax = {}'.format(ymin, ymax))


fig.subplots_adjust(right=0.98,left=0.05,bottom=0.15,top = 0.98,wspace=0.3)
# fig2.subplots_adjust(right=0.98,left=0.2,bottom=0.2,top = 0.98)
# outfn = home+'plots/Feb2023_moorings_fit.png'
# plt.savefig(outfn)
# print(f'Saved to {outfn}')
plt.show(block=False)
plt.pause(0.1)
# plt.close()


