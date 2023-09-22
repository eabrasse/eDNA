"""
Plot results of a particle tracking experiment.
do I want to plot all 40 in one figure? or plot them in 40 figures? or a subset in one figure?
"""


import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import matplotlib
import pickle
import pytz
from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LogNorm

Ldir = Lfun.Lstart()

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

#custom sequential colormap
tab10 = plt.get_cmap('tab10',10)
mycolor = [p for p in tab10(0)]
colors = [(mycolor[0], mycolor[1], mycolor[2],c) for c in np.linspace(0,1,100)]
newBlues = matplotlib.colors.LinearSegmentedColormap.from_list('newBlue', colors, N=100)

#custom sequential colormap
mycolor = [p for p in tab10(1)]
colors = [(mycolor[0], mycolor[1], mycolor[2],c) for c in np.linspace(0,1,100)]
newOranges = matplotlib.colors.LinearSegmentedColormap.from_list('newOranges', colors, N=100)

# #custom diverging colormap
# ncolors = 256
# ncolors_half = int(ncolors/2)
# myblue = [p for p in matplotlib.colors.to_rgba('blue')]
# myred = [p for p in matplotlib.colors.to_rgba('red')]
# blue_array = [(myblue[0], myblue[1], myblue[2],c) for c in np.linspace(0,1,ncolors_half)[::-1]]
# red_array = [(myred[0], myred[1], myred[2],c) for c in np.linspace(0,1,ncolors_half)]
# color_array = np.concatenate((blue_array,red_array),axis=0)
# newbwr = matplotlib.colors.LinearSegmentedColormap.from_list('newbwr', color_array, N=100)

# # Choose an experiment and release to plot.
in_dir0 = Ldir['LOo'] / 'tracks'
data_fn = in_dir0 / 'all_3d_hc_dolph_hist_by_age.p'

D = pickle.load(open(data_fn,'rb'))

# gather some fields, for convenience
lonp = D['lon'][:]
latp = D['lat'][:]
hh = D['h'][:]
maskr = D['mask'][:]

hr_list = [1, 12, 24, 36, 48]
# hr_list = [hr*3600 for hr in hr_list] #translate to seconds

# mask path
# mask_fn = '/data2/pmr4/eab32/LO_output/extract/hc11_v1_uu0k/wetdry_mask_tides.p'
# Dmask = pickle.load(open(mask_fn,'rb'))
# wd_mask = Dmask['wetdry_mask_rho'][:]
# wd_dt_list = Dmask['dt_list'][:]
# wd_dt_list = [pytz.utc.localize(wd_dt) for wd_dt in wd_dt_list]

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

# MAP
# set domain limits
pad = .01
# plot region around delta pier, using release point as indicator
aa = [lon0-pad, lon0+pad,
   lat0-pad, lat0+pad]
 
#identify grid edge limits for making mask      
AA = [lonp[0,0], lonp[0,-1],
        latp[0,0], latp[-1,0]]
xx, yy = np.meshgrid(D['bin_lon_edges'][:-1]+0.5*(D['bin_lon_edges'][1]-D['bin_lon_edges'][0]),
    D['bin_lat_edges'][:-1]+0.5*(D['bin_lat_edges'][1]-D['bin_lat_edges'][0]))
# vmax = 75

# PLOTTING
plt.close('all')
pfun.start_plot(figsize=(10,8))
fig,axes = plt.subplots(2,len(hr_list))
count = 0
for hr in hr_list:
    
    ax_mean = axes[0,count]
    ax_var = axes[1,count]
    
    zm = -np.ma.masked_where(maskr==0, hh)
    
    for ax in ax_mean, ax_var:
        # mask land
        ax.pcolormesh(lonp,latp,maskr,cmap=cmap_mask,shading='nearest',zorder=0)
        # indicate release point
        ax.plot(lon0, lat0, marker='*',mec='k',mfc='yellow', markersize=10,alpha=1,zorder=500)
    
        ax.axis(aa)
        pfun.dar(ax)
    titlestr = f'{hr:} hours'
    
    particle_mean = np.sum(D['hist'][:,hr,:,:],axis=0) #sum all releases for 'hr' hours since release
    particle_var = np.std(D['hist'][:,hr,:,:],axis=0)**2 #estimate variance across releases for 'hr' hours since release
    
    pm0=ax_mean.pcolormesh(xx,yy,particle_mean,norm=matplotlib.colors.LogNorm(vmin=1e1,vmax=1e8),cmap=newBlues)
    pv0=ax_var.pcolormesh(xx,yy,particle_var,norm=matplotlib.colors.LogNorm(vmin=1e1,vmax=1e16),cmap=newOranges)
    
    if hr==1:
        titlestr = titlestr[:-1]
        pv = pv0
        pm = pm0
    ax_mean.text(0.9,0.9,titlestr,transform=ax_mean.transAxes,ha='right',va='top')
    
    count+=1
    
axes[0,0].text(0.1,0.9,'MEAN',transform=axes[0,0].transAxes)
cbaxes = inset_axes(axes[0,-1], width="4%", height="60%", loc='center right',bbox_transform=axes[0,-1].transAxes,bbox_to_anchor=(0.15,0.,1,1))
cb = fig.colorbar(pm, cax=cbaxes, orientation='vertical')
cb.set_label('mean particles')

axes[1,0].text(0.1,0.9,'VARIANCE',transform=axes[1,0].transAxes)
cbaxes = inset_axes(axes[1,-1], width="4%", height="60%", loc='center right',bbox_transform=axes[1,-1].transAxes,bbox_to_anchor=(0.15,0.,1,1))
cb = fig.colorbar(pv, cax=cbaxes, orientation='vertical')
cb.set_label('variance particles')

for ax in axes[:,0]:
    ax.set_ylabel('Latitude')
for ax in axes[:,1:].flatten():
    ytl = ['' for yt in ax.get_yticklabels()]
    ax.set_yticklabels(ytl)

for ax in axes[1,:]:
    ax.set_xlabel('\n Longitude')
for ax in axes[0,:]:
    xtl = ['' for xt in ax.get_xticklabels()]
    ax.set_xticklabels(xtl)

# plt.tight_layout()
outfn = '/data2/pmr4/eab32/etools/plots/hc_dolph_3d_hist_by_age.png'
plt.savefig(outfn)
pfun.end_plot()

print('saved to {}'.format(outfn))
    
    

