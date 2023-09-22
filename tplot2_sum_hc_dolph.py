"""
Plot results of a particle tracking experiment.
do I want to plot all 40 in one figure? or plot them in 40 figures? or a subset in one figure?
"""


import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pickle
import pytz
from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

#custom sequential colormap
paired = plt.get_cmap('Paired',12)
mygreen = [p for p in paired(3)]
colors = [(mygreen[0], mygreen[1], mygreen[2],c) for c in np.linspace(0,1,100)]
newBlues = matplotlib.colors.LinearSegmentedColormap.from_list('mycmap', colors, N=100)

#custom diverging colormap
ncolors = 256
ncolors_half = int(ncolors/2)
myblue = [p for p in matplotlib.colors.to_rgba('blue')]
myred = [p for p in matplotlib.colors.to_rgba('red')]
blue_array = [(myblue[0], myblue[1], myblue[2],c) for c in np.linspace(0,1,ncolors_half)[::-1]]
red_array = [(myred[0], myred[1], myred[2],c) for c in np.linspace(0,1,ncolors_half)]
color_array = np.concatenate((blue_array,red_array),axis=0)
newbwr = matplotlib.colors.LinearSegmentedColormap.from_list('newbwr', color_array, N=100)

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

# MAP
# set domain limits
pad = .01
# plot region around delta pier, using release point as indicator
aa = [lon0-pad, lon0+pad,
   lat0-pad, lat0+pad]
 
#identify grid edge limits for making mask      
AA = [lonp[0,0], lonp[0,-1],
        latp[0,0], latp[-1,0]]

nbins = 50
bin_lon_edges=np.linspace(aa[0], aa[1],nbins+1)
bin_lat_edges=np.linspace(aa[2], aa[3],nbins+1)
xx, yy = np.meshgrid(bin_lon_edges[:-1]+0.5*(bin_lon_edges[1]-bin_lon_edges[0]),bin_lat_edges[:-1]+0.5*(bin_lat_edges[1]-bin_lat_edges[0]))

vmax = 75
const_denom=0
tvary_denom = 0
pulse = {'desc':'Pulse release','hist':0}
const = {'desc':'Hourly releases\nNo decay\nConst shedding','hist':0}
decay = {'desc':'Hourly releases\nDecay\nConst shedding','hist':0}
tvary = {'desc':'Hourly releases\nDecay\nVariable shedding.','hist':0}
NP = 100000
for release in release_list:
    rel = D[release]

    lon = rel['lon'][:]
    lat = rel['lat'][:]

    # make a mask that is False from the time a particle first leaves the domain
    # and onwards

    ib_mask = np.ones(lon.shape, dtype=bool)
    ib_mask[lon < AA[0]] = False
    ib_mask[lon > AA[1]] = False
    ib_mask[lat < AA[2]] = False
    ib_mask[lat > AA[3]] = False

    # and apply the mask to lon and lat
    lon[~ib_mask] = np.nan
    lat[~ib_mask] = np.nan
    
    hist = np.histogram2d(lon[:],lat[:],bins=[bin_lon_edges,bin_lat_edges])
    
    # add to hist piles
    if release=='rel_2023.01.31_sh1':
        pulse['hist'] = hist[0].T
    const['hist'] += hist[0].T
    decay['hist'] += rel['decay']*hist[0].T
    tvary['hist'] += rel['decay']*rel['C0']*hist[0].T
    
    const_denom += NP
    tvary_denom += rel['C0']*NP

#normalize
pulse['hist'] = pulse['hist']/NP
const['hist'] = const['hist']/const_denom
decay['hist'] = decay['hist']/const_denom
tvary['hist'] = tvary['hist']/tvary_denom

release_type_list = [tvary,decay,const,pulse]
Nrt = len(release_type_list)

# PLOTTING
plt.close('all')
pfun.start_plot(figsize=(14,8))
fig,axes = plt.subplots(2,Nrt)

vmin = 0
vmax_list = [np.percentile(rt['hist'][rt['hist']>0],90) for rt in release_type_list]
vmax = np.max(vmax_list)
vmax2 = 0.001
count=0
for release_type in release_type_list:
    ax0 = axes[0,count]
    ax1 = axes[1,count]
    
    for ax in [ax0,ax1]:
        #add mask and wetdry mask
        ax.pcolormesh(lonp,latp,maskr,cmap=cmap_mask,shading='nearest',zorder=0)
        ax.contour(lonp,latp,wd_mask[wd_ind,:,:],levels=[0.5],colors=['k'],linewidths=[1.5])

        ax.axis(aa)
        pfun.dar(ax)
        ax.plot(lon0, lat0, marker='*',mec='k',mfc='yellow', markersize=10,alpha=1,zorder=500)
    
    p=ax0.pcolormesh(xx,yy,release_type['hist'],vmax=vmax,vmin=vmin,cmap=newBlues)
    p2 = ax1.pcolormesh(xx,yy,release_type['hist']-tvary['hist'],vmax=vmax2,vmin=-vmax2,cmap=newbwr)


    if count==0:
        ax0.set_ylabel('Latitude')
    else:
        ax0.set_yticklabels([''])
    ax0.set_title(release_type['desc'])
    ax0.set_xticklabels([''])
    ax1.set_xlabel('Longitude')
    
    count+=1
# plt.colorbar(p)

cbaxes = inset_axes(axes[0,-1], width="4%", height="60%", loc='center right',bbox_transform=axes[0,-1].transAxes,bbox_to_anchor=(0.15,0.,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
cb.set_label('fraction of released particles')

cbaxes = inset_axes(ax1, width="4%", height="60%", loc='center right',bbox_transform=ax1.transAxes,bbox_to_anchor=(0.15,0.,1,1))
cb = fig.colorbar(p2, cax=cbaxes, orientation='vertical')
cb.set_label('diff from variable shedding')

#plt.show()
outfn = '/data2/pmr4/eab32/etools/plots/hc_dolph_3d_compare_release_typesB.png'
plt.savefig(outfn)
pfun.end_plot()

print('saved to {}'.format(outfn))

