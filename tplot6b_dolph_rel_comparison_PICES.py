"""
Plot results of a particle tracking experiment.
do I want to plot all 40 in one figure? or plot them in 40 figures? or a subset in one figure?
"""

import os
import sys
# pth = os.path.abspath('/data2/pmr4/eab32/etools/')
# if pth not in sys.path:
#     sys.path.append(pth)
import efun
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pickle
import pytz
from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
import string
import cmocean as cmo

atoz = string.ascii_lowercase
Ldir = Lfun.Lstart()

tab10 = plt.get_cmap('tab10',10)

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

# #custom sequential colormap
# paired = plt.get_cmap('Paired',12)
# mygreen = [p for p in paired(3)]
# colors = [(mygreen[0], mygreen[1], mygreen[2],c) for c in np.linspace(0,1,100)]
# newBlues = matplotlib.colors.LinearSegmentedColormap.from_list('mycmap', colors, N=100)


# # Choose an experiment and release to plot.
in_dir0 = Ldir['LOo'] / 'tracks'
data_fn = in_dir0 / 'all_3d_hc_dolph_releases.p'

D = pickle.load(open(data_fn,'rb'))

# gather some fields, for convenience
lonp = D['metadata']['lon'][:]
latp = D['metadata']['lat'][:]
hh = D['metadata']['h'][:]
maskr = D['metadata']['mask'][:]


key_list = D.keys()
release_list = [key for key in key_list if key[:3]=='rel']


# mask path
mask_fn = '/data2/pmr4/eab32/LO_output/extract/hc11_v1_uu0k/wetdry_mask_tides.p'
Dmask = pickle.load(open(mask_fn,'rb'))
wd_mask = Dmask['wetdry_mask_rho'][:]
wd_dt_list = Dmask['dt_list'][:]
wd_dt_list = [pytz.utc.localize(wd_dt) for wd_dt in wd_dt_list]
wd_td_list = [np.abs(wd_dt - D['metadata']['refT']) for wd_dt in wd_dt_list]
wd_ind = wd_td_list.index(min(wd_td_list))

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

# # MAP
# # set domain limits
# pad = .01
# # plot region around delta pier, using release point as indicator
# aa = [lon0-pad, lon0+pad,
#    lat0-pad, lat0+pad]
#
# #identify grid edge limits for making mask
# AA = [lonp[0,0], lonp[0,-1],
#         latp[0,0], latp[-1,0]]
# zref = -2
#
# nbins = 100
# bin_lon_edges=np.linspace(aa[0], aa[1],nbins+1)
# bin_lat_edges=np.linspace(aa[2], aa[3],nbins+1)
# xx, yy = np.meshgrid(bin_lon_edges[:-1]+0.5*(bin_lon_edges[1]-bin_lon_edges[0]),bin_lat_edges[:-1]+0.5*(bin_lat_edges[1]-bin_lat_edges[0]))
# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

# MAP
xgrid,ygrid = efun.ll2xy(lonp,latp,lon0,lat0)
# set domain limits
pad = .01
# plot region around delta pier, using release point as indicator
aa = [lon0-pad, lon0+pad,
   lat0-pad, lat0+pad]

aax,aay =  efun.ll2xy(np.array(aa[:2]),np.array(aa[2:]),lon0,lat0)
aaxy = [aax[0],aax[1],aay[0],aay[1]]

#identify grid edge limits for making mask      
AA = [lonp[0,0], lonp[0,-1],
        latp[0,0], latp[-1,0]]
zref = -2

nbins = 100
bin_x_edges=np.linspace(aax[0], aax[1],nbins+1)
bin_y_edges=np.linspace(aay[0], aay[1],nbins+1)
xx, yy = np.meshgrid(bin_x_edges[:-1]+0.5*(bin_x_edges[1]-bin_x_edges[0]),bin_y_edges[:-1]+0.5*(bin_y_edges[1]-bin_y_edges[0]))
# dxbin = bin_x_edges[1]-bin_x_edges[0]
# dybin = bin_y_edges[1]-bin_y_edges[0]
# dzbin = np.abs(zref)

vmax = 75
const_denom=0
tvary_denom = 0
pulse = {'desc':'Pulse','hist_nodecay':0,'hist_decay':0,'denom':0}
const = {'desc':'Constant','hist_nodecay':0,'hist_decay':0,'denom':0}
tvary = {'desc':'Time-varying','hist_nodecay':0,'hist_decay':0,'denom':0}
# NP = 100000
for release in release_list:
    rel = D[release]
    # print('deltaT={}'.format(rel['deltaT']))
    lon = rel['lon'][:]
    lat = rel['lat'][:]
    z = rel['z'][:]

    # make a mask that is False from the time a particle first leaves the domain
    # and onwards

    ib_mask = np.ones(lon.shape, dtype=bool)
    ib_mask[lon < AA[0]] = False
    ib_mask[lon > AA[1]] = False
    ib_mask[lat < AA[2]] = False
    ib_mask[lat > AA[3]] = False
    ib_mask[z<zref] = False

    # and apply the mask to lon and lat
    lon[~ib_mask] = np.nan
    lat[~ib_mask] = np.nan
    
    # hist = np.histogram2d(lon[:],lat[:],bins=[bin_lon_edges,bin_lat_edges])
    xp,yp = efun.ll2xy(lon,lat,lon0,lat0)
    
    hist = np.histogram2d(xp,yp,bins=[bin_x_edges,bin_y_edges])
    
    # add to hist piles
    if release==release_list[0]:
        NP = lon.shape[0]
    if rel['deltaT']==(20*3600):
        pulse['hist_nodecay'] = hist[0].T
        pulse['hist_decay'] = rel['decay']*hist[0].T
        pulse['denom']+=NP
    const['hist_nodecay'] += hist[0].T
    const['hist_decay'] += rel['decay']*hist[0].T
    tvary['hist_nodecay'] += rel['C0']*hist[0].T
    tvary['hist_decay'] += rel['decay']*rel['C0']*hist[0].T
    
    const['denom'] += NP
    tvary['denom'] += rel['C0']*NP

#normalize
release_type_list = [tvary,const,pulse]
decay_key_list = ['hist_nodecay','hist_decay']
decay_label_list = ['','+ Decay']
for reltype in release_type_list:
    for key in decay_key_list:
        reltype[key] = reltype[key]/reltype['denom']

Nrt = len(release_type_list)
Ndk = len(decay_key_list)
vmin = 0
vmax_list = [np.percentile(100*rt['hist_nodecay'][rt['hist_nodecay']>0],90) for rt in release_type_list]
vmax = np.max(vmax_list)
vmax2 = 0.001
# xymax = 100*max([rt['hist_nodecay'].max() for rt in release_type_list])
xymax = 100*max([np.percentile(rt['hist_nodecay'],98) for rt in release_type_list])
xmean = tvary['hist_decay'].mean()
# sort indexes
tv_inds = tvary['hist_decay'].flatten().argsort()
xsort = 100*tvary['hist_decay'].flatten()[tv_inds]

# PLOTTING
plt.close('all')
fw, fh = efun.gen_plot_props()
figsize = (fw*3,fh)
fig1,axes1 = plt.subplots(1,Nrt,figsize=figsize)
# fig2,axes2 = plt.subplots(1,Nrt,figsize=figsize)
count=0
for rt in range(Nrt):
    release_type = release_type_list[rt]
    # for dk in range(Ndk):
    dk = 1
    key = decay_key_list[dk]
    
    # map of distribution
    ax1 = axes1[rt]
    
    # add context - shoreline, etc
    ax1.pcolormesh(xgrid,ygrid,maskr,cmap=cmap_mask,shading='nearest',zorder=0)
    ax1.contour(xgrid,ygrid,wd_mask[wd_ind,:,:],levels=[0.5],colors=['k'],linewidths=[1.5])
    ax1.axis(aaxy)
    # pfun.dar(ax1)
    ax1.plot([0], [0], marker='*',mec='k',mfc='yellow', markersize=10,alpha=1,zorder=500)
    # plot particle heatmap
    release_type[key][release_type[key]==0]=np.nan
    p=ax1.pcolormesh(xx,yy,100*release_type[key],vmax=vmax,vmin=vmin,cmap=cmo.cm.matter)
    # add axis labels
    if rt>0:
        ax1.set_yticklabels([''])
    else:
        ax1.set_ylabel('Dist from pen (m)')
    ax1.set_xlabel('Dist from pen (m)')
    plt.setp( ax1.xaxis.get_majorticklabels(), rotation=15, ha="right",rotation_mode='anchor')
    ax1.ticklabel_format(axis='x',style='plain',useOffset=False)
    # ax1.set_title(release_type['desc'])
    

    # # now plot comparison of each option with Time-Varying decay
    # ax2 = axes2[rt]
    # ax2.scatter(100*tvary['hist_decay'].flatten(),100*release_type[key].flatten(),color=mygreen,s=10,alpha=0.35,zorder=30)
    # ax2.set_aspect(1)
    # # xx = tvary['hist_decay'].flatten().sort()
    # ax2.plot([0,xymax],[0,xymax],linestyle='solid',color='gray',zorder=25)
    # ax2.axis([0,xymax,0,xymax])
    # rmse = np.sqrt((np.sum((tvary['hist_decay'].flatten()-release_type[key].flatten())**2))/len(tvary['hist_decay'].flatten()))
    # nrmse = rmse/xmean
    #
    # ax2.set_ylabel(release_type['desc']+decay_label_list[dk])
    # ax2.set_xlabel(tvary['desc']+decay_label_list[1])
    # ysort = [100*release_type[key].flatten()[tv_ind] for tv_ind in tv_inds]
    # pfit = np.polyfit(xsort,ysort,1)
    # xfit = np.linspace(0,xymax,10)
    # yfit = pfit[0]*xfit +pfit[1]
    # ax2.plot(xfit,yfit,linestyle='dashed',color='magenta',zorder=26)
    # R = np.corrcoef(xsort,ysort)[0,1]
    # Rsquared = R**2
    # wss = efun.willmott(xsort,ysort)
    #
    # ax2.text(0.98,0.1,
    #     f'RMSE = {rmse:0.2}\nNRSME = {nrmse:0.2}\nlinear fit = {pfit[0]:0.2}x+{pfit[1]:0.2}\nR2={Rsquared:0.2}\nWSS={wss:0.2}',
    #     transform=ax2.transAxes,ha='right',va='bottom',color='k',zorder=50,fontsize=10)
    
    ax1.set_aspect(1)
    ax1.text(0.1,0.9,'{:}) {:}'.format(atoz[count],release_type['desc']),transform=ax1.transAxes,ha='left',va='top',color='black')
    # ax2.text(0.1,0.9,'{:}) {:}{:}'.format(atoz[count],release_type['desc'],decay_label_list[dk]),transform=ax2.transAxes,ha='left',va='top',color='black')
    
    count+=1



cbaxes = inset_axes(axes1[-1], width="4%", height="60%", loc='center right',bbox_transform=axes1[-1].transAxes,bbox_to_anchor=(0.15,0.,1,1))
cb = fig1.colorbar(p, cax=cbaxes, orientation='vertical')
cb.set_label(r'$\%$ released particles')

#plt.show()
fig1.subplots_adjust(bottom=0.2,top=0.98,left=0.08,right=0.88,wspace=0.05)
outfn1 = '/data2/pmr4/eab32/etools/plots/hc_dolph_3d_compare_releases_heatmapC.png'
fig1.savefig(outfn1)
print('saved to {}'.format(outfn1))
# fig2.subplots_adjust(bottom=0.08,top=0.98,left=0.08,right=0.98,hspace=0.3)
# outfn2 = '/data2/pmr4/eab32/etools/plots/hc_dolph_3d_compare_releases_scatterB.png'
# fig2.savefig(outfn2)
# print('saved to {}'.format(outfn2))


plt.close('all')
