"""
Plot results of a particle tracking experiment.
do I want to plot all 40 in one figure? or plot them in 40 figures? or a subset in one figure?
"""

import os
import sys
pth = os.path.abspath('/data2/pmr4/eab32/etools/')
if pth not in sys.path:
    sys.path.append(pth)
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

atoz = string.ascii_lowercase
Ldir = Lfun.Lstart()

tab10 = plt.get_cmap('tab10',10)

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

#custom sequential colormap
paired = plt.get_cmap('Paired',12)
mygreen = [p for p in paired(3)]
colors = [(mygreen[0], mygreen[1], mygreen[2],c) for c in np.linspace(0.25,1,99)]
colors.insert(0,(mygreen[0], mygreen[1], mygreen[2],0))
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
dxbin = bin_x_edges[1]-bin_x_edges[0]
dybin = bin_y_edges[1]-bin_y_edges[0] 
dzbin = np.abs(zref)

vmax = 75

particle_bin = 0
# NP = 100000
for release in release_list:
    rel = D[release]

    lon = rel['lon'][:]
    lat = rel['lat'][:]
    z = rel['z'][:]
    
    if rel['deltaT']==0:
        # get number of particles
        particle_rel = lon.shape[0]
        # derive the area over which they were released
        xp0,yp0 = efun.ll2xy(lon,lat,lon0,lat0)
        dxrel = np.abs(xp0.max()-xp0.min())
        dyrel = np.abs(yp0.max()-yp0.min())
        dzrel = 1 #all released at surface, soooo....

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
    
    xp,yp = efun.ll2xy(lon,lat,lon0,lat0)
    
    hist = np.histogram2d(xp,yp,bins=[bin_x_edges,bin_y_edges])
    
    # add to hist piles
    # particle_bin += rel['decay']*hist[0].T
    particle_bin += hist[0].T

#normalize
vol_rel = dxrel*dyrel*dzrel
particle_conc_rel = particle_rel/vol_rel

vol_bin = dxbin*dybin*dzbin
particle_conc_bin = particle_bin/vol_bin

DNA_conc_rel_list = [1e2,1e3,1e4,1e5]

nDNA_conc_rel = len(DNA_conc_rel_list)
plt.close('all')
fw,fh=efun.gen_plot_props()
figsize = (fw*4,fh*1.5)
fig,axs = plt.subplots(1,nDNA_conc_rel,figsize=figsize)
vmin = 0
# vmax = max(DNA_conc_rel_list) * np.percentile(particle_conc_bin,90)/particle_conc_rel
vmax = 5

for i in range(nDNA_conc_rel):
    DNA_conc_rel = DNA_conc_rel_list[i]
    DNA_conc_bin = DNA_conc_rel * particle_conc_bin/particle_conc_rel
    ax = axs[i]
    ax.pcolormesh(xgrid,ygrid,maskr,cmap=cmap_mask,shading='nearest',zorder=0)
    ax.contour(xgrid,ygrid,wd_mask[wd_ind,:,:],levels=[0.5],colors=['k'],linewidths=[1.5])
    ax.axis(aaxy)
    ax.plot([0],[0],marker='*',mec='k',mfc='yellow',markersize=15,alpha=1,zorder=500)
    p=ax.pcolormesh(xx,yy,DNA_conc_bin,vmax=vmax,vmin=vmin,cmap=newBlues)
    if DNA_conc_bin.max()>5:
        ax.contour(xx,yy,DNA_conc_bin,levels=[5],colors=['m'],linewidths=[1],linestyles=['solid'],zorder=1000)
    ax.set_xlabel('Dist from pen (m)')
    ax.set_ylabel('Dist from pen (m)')
    ax.text(0.1,0.95,'DNA conc near source\n'+f'{int(DNA_conc_rel):}'+r' copies $\mu\mathrm{L}^{-1}$',transform=ax.transAxes,ha='left',va='top')

cbaxes = inset_axes(axs[-1], width="4%", height="60%", loc='center right',bbox_transform=axs[-1].transAxes,bbox_to_anchor=(0.15,0.,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
cb.set_label(r'copies $\mu\mathrm{L}^{-1}$')

#plt.show()
plt.subplots_adjust(left=0.1,right=0.9,wspace=0.4)
outfn = '/data2/pmr4/eab32/etools/plots/hc_dolph_LOD.png'
fig.savefig(outfn)
print('saved to {}'.format(outfn))

plt.close('all')

#now plot area vs release count
npoints = 100
DNA_conc_rel_lgsp = np.logspace(2,5,npoints)
LOD_area = np.zeros((npoints))
for i in range(npoints):
    DNA_conc_bin = DNA_conc_rel_lgsp[i] * particle_conc_bin/particle_conc_rel
    count_bins_above_LOD = np.sum(DNA_conc_bin>5)
    LOD_area[i] = dxbin*dybin*count_bins_above_LOD
fig= plt.figure(figsize=(1.5*fw,1.25*fh))
ax=fig.gca()
ax.plot(DNA_conc_rel_lgsp,LOD_area,color=tab10(0),lw=3)
ax.set_xscale('log')
ax.set_xlabel(r'DNA concentration near source (copies $\mu\mathrm{L}^{-1}$)')
ax.set_ylabel(r'Area above LOD ($\mathrm{m}^{2}$)')
ax.grid()

plt.subplots_adjust(left=0.23,right=0.95,top=0.98,bottom=0.18)
outfn = '/data2/pmr4/eab32/etools/plots/hc_dolph_LOD_area.png'
fig.savefig(outfn)
print('saved to {}'.format(outfn))
plt.close('all')