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
from scipy.stats import linregress
import cmocean as cmo
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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
moor_fn = home+'data/Feb2023_DNA_moorings_extended_const_only_noscale.p'

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

td_list = [np.abs(dt-sample_dt) for dt in dt_list]
sample_dt_ind = np.argmin(td_list)


fw,fh = efun.gen_plot_props()
fig2 = plt.figure(figsize=(fw,fh))
ax_samp = fig2.gca()

# # get grid from random file
grid_fn= home+ 'data/all_3d_hc_dolph_releases.p'

D = pickle.load(open(grid_fn,'rb'))

# gather some fields, for convenience
lonr = D['metadata']['lon'][:]
latr = D['metadata']['lat'][:]
maskr = D['metadata']['mask'][:]


props = dict(boxstyle='round', facecolor='white',edgecolor='none',alpha=0.8)

moor_count=0
mooring_axes = []
particles = np.zeros((nmoor-1))
samples = np.zeros((nmoor-1))
for ind in ns_inds:

    moor = moor_list[ind]
    if moor['label'] in mll2ref.keys():
        reference=mll2ref[moor['label']]
        gg = gq[gq['reference']==reference]
        sampledata = gg['PB_quantity'].values
        mean_sampledata = np.nanmean(sampledata)
        
        particles[moor_count] = moor['const_DNA_bin'][sample_dt_ind]
        nsamples = len(sampledata)
        samples[moor_count] = mean_sampledata
        
        ax_samp.plot(particles[moor_count],samples[moor_count],marker='d',mfc=rainbow(moor_count),mec='k',linestyle='None',markersize=10,alpha=1.0,zorder=15)
        ax_samp.plot(particles[moor_count]*np.ones((nsamples)),sampledata,marker='d',mfc=rainbow(moor_count),mec='None',linestyle='None',markersize=6,alpha=0.5,zorder=5)
        moor_count+=1
        
pinds = np.argsort(particles)
x0 = particles[pinds]
y0 = samples[pinds]
# p = np.polyfit(x0,y0,1)
slope, intercept, r_value, p_value, std_err = linregress(x0,y0)
xx = np.linspace(x0.min(),x0.max(),100)
yy = xx*slope+intercept
ax_samp.plot(xx,yy,linestyle='dashed',color='k',zorder=10)
# ax_samp.text(0.5,0.5,'best fit line = {:0.0f}x+{:0.1f}'.format(slope,intercept),transform=ax_samp.transAxes,ha='center',va='center')
ax_samp.set_xlabel('scaled particles')
ax_samp.set_ylabel('sampled DNA concentration')
ax_samp.set_xscale('log')
ax_samp.set_yscale('log')
ax_samp.text(0.4,0.65,'Linear fit',transform=ax_samp.transAxes,rotation=30)
ax_samp.text(0.4,0.4,r'$\mathrm{R}^{2}$ '+'= {:.4}'.format(r_value**2)+'\n'+'p={:.2E}'.format(p_value),transform=ax_samp.transAxes)

# MAP
fig = plt.figure(figsize=(14, 8))
gs = GridSpec(nmoor-1+1,3)

ax = plt.subplot(gs[:,0])
xgrid,ygrid = efun.ll2xy(lonr,latr,lon0,lat0)
#~~~~~~~~~~~~
#reprint code to calculate bins because I forgot to output them
# set domain limits
pad = .01
# plot region around delta pier, using release point as indicator
aa = [lon0-pad, lon0+pad,
   lat0-pad, lat0+pad]
aax,aay =  efun.ll2xy(np.array(aa[:2]),np.array(aa[2:]),lon0,lat0)
aaxy = [aax[0],aax[1],aay[0],aay[1]]
nbins = 100
bin_x_edges=np.linspace(aax[0], aax[1],nbins+1)
xbins = 0.5*(bin_x_edges[:-1]+bin_x_edges[1:])
bin_y_edges=np.linspace(aay[0], aay[1],nbins+1)
ybins = 0.5*(bin_y_edges[:-1]+bin_y_edges[1:])
#~~~~~~~~~~~~~
ax.pcolormesh(xgrid,ygrid,maskr,cmap=cmap_mask,shading='nearest',zorder=100)
p=ax.pcolormesh(xbins,ybins,slope*moor_dict['const_DNA_bin'][sample_dt_ind,:].T,cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(),zorder=50,shading='nearest')
ax.contour(xgrid,ygrid,maskr,levels=[0.5],colors=['k'],linewidths=[1.5],zorder=150)
ax.set_xlabel('Dist from pen (m)')
ax.set_ylabel('Dist from pen (m)')
cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower right',bbox_transform=ax.transAxes,bbox_to_anchor=(-0.2,0.,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
cb.set_label('scaled particles')
moor_count=0
for ind in ns_inds:

    moor = moor_list[ind]
    if moor['label'] in mll2ref.keys():
        
    
        x,y = efun.ll2xy(moor['lon'],moor['lat'],lon0,lat0)
        ax.plot(x,y,marker='d',mec='k',mfc = rainbow(moor_count),markersize=15,zorder=200)

        axm = plt.subplot(gs[moor_count,1:])

        axm.plot(dt_list,slope*moor['const_DNA_bin'],linestyle='dashed',color=rainbow(moor_count))


    
        axm.text(0.05,0.8,'{}) {}'.format(atoz[moor_count],moor['label']),color='k',transform=axm.transAxes,ha='left',va='top')
        axm.grid()

        # get universal ylims
        if moor_count==0:
            ymax = np.nanmax(slope*moor['const_DNA_bin'])
            ymin = np.nanmin(slope*moor['const_DNA_bin'][moor['const_DNA_bin']>0])
        ymax = max([np.nanmax(slope*moor['const_DNA_bin']),ymax])
        ymin = min([np.nanmin(slope*moor['const_DNA_bin'][moor['const_DNA_bin']>0]),ymin])

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
        # if moor_count==0:
            # axm.text(0.8,0.95,'Solid = time-varying\nDashed = constant\nDots = samples\nDotted line = sample mean',transform=axm.transAxes,ha='right',va='top',bbox=props,fontsize=10)

        moor_count+=1


axp = plt.subplot(gs[-1,1:])
axp.fill_between(dt_list,np.zeros((len(dt_list))),moor_dict['active_particles'],color='k',alpha=0.5)
axp.set_xlabel('Date & time (UTC)')
axp.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d %H:%m"))
plt.setp( axp.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
axp.text(0.05,0.8,'{}) {}'.format(atoz[moor_count],'Active particles'),color='k',transform=axp.transAxes,ha='left',va='top')
axp.grid()
dom_part = np.sum(np.sum(moor_dict['const_DNA_bin'][:],axis=2),axis=1)
axp.fill_between(dt_list,dom_part,color='k',alpha=0.5)

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
fig2.subplots_adjust(right=0.98,left=0.2,bottom=0.2,top = 0.98)
# outfn = home+'plots/Feb2023_moorings_obs.png'
# plt.savefig(outfn)
# print(f'Saved to {outfn}')
plt.show(block=False)
plt.pause(0.1)
# plt.close()


