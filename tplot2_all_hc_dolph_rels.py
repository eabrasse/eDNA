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
Ldir = Lfun.Lstart()

# useful plotting tools, colormap and test box props
landcol = 'lightgray'
seacol = 'white'

cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])
paired = plt.get_cmap('Paired',12)
myblue = [p for p in paired(1)]
colors = [(myblue[0], myblue[1], myblue[2],c) for c in np.linspace(0,1,100)]
newBlues = matplotlib.colors.LinearSegmentedColormap.from_list('mycmap', colors, N=100)

# # Choose an experiment and release to plot.
in_dir0 = Ldir['LOo'] / 'tracks'
data_fn = in_dir0 / 'all_3d_hc_dolph_releases.p'
# exp_name = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
#     itext='** Choose experiment from list **', last=False)
# rel = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
#     itext='** Choose item from list **', last=False)

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

vmax = 75
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
    # NTS, NPS = lon.shape
    # for pp in range(NPS):
    #     tt = np.argwhere(ib_mask[:,pp]==False)
    #     if len(tt) > 0:
    #         ib_mask[tt[0][0]:, pp] = False

    # and apply the mask to lon and lat
    lon[~ib_mask] = np.nan
    lat[~ib_mask] = np.nan

    # PLOTTING
    plt.close('all')
    pfun.start_plot(figsize=(10,8))
    fig = plt.figure()


    
    # ax = fig.add_subplot(121)
    ax = fig.gca()
    zm = -np.ma.masked_where(maskr==0, hh)
    # plt.pcolormesh(lonp, latp, zm, vmin=-100, vmax=0,
        # cmap='terrain', alpha=.25)
    # pfun.add_coast(ax)
    #add mask and wetdry mask
    wd_td_list = [np.abs(wd_dt - D['metadata']['refT']) for wd_dt in wd_dt_list]
    wd_ind = wd_td_list.index(min(wd_td_list))
    ax.pcolormesh(lonp,latp,maskr,cmap=cmap_mask,shading='nearest',zorder=0)
    ax.contour(lonp,latp,wd_mask[wd_ind,:,:],levels=[0.5],colors=['k'],linewidths=[1.5])
    
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(rel['T0'])
    # add the tracks (packed [time, particle])
    # regular spaghetti plots
    # ax.plot(lon, lat, '-k', linewidth=.2)
    ax.plot(lon0, lat0, marker='*',mec='k',mfc='yellow', markersize=10,alpha=1,zorder=500)
    # ax.plot(lon[-1,:], lat[-1,:], 'or', alpha=.3)
    nbins = 50
    bin_lon_edges=np.linspace(aa[0], aa[1],nbins+1)
    bin_lat_edges=np.linspace(aa[2], aa[3],nbins+1)
    xx, yy = np.meshgrid(bin_lon_edges[:-1]+0.5*(bin_lon_edges[1]-bin_lon_edges[0]),bin_lat_edges[:-1]+0.5*(bin_lat_edges[1]-bin_lat_edges[0]))
    hist = np.histogram2d(lon[:],lat[:],bins=[bin_lon_edges,bin_lat_edges])
    # vmax=np.percentile(hist[0][hist[0]>0],90)
    p=ax.pcolormesh(xx,yy,hist[0].transpose(),vmax=vmax,cmap=newBlues)
    plt.colorbar(p)


    #plt.show()
    outfn = '/data2/pmr4/eab32/etools/plots/hc_dolph_3d_heatmaps/'+release+'.png'
    plt.savefig(outfn)
    pfun.end_plot()

    print('saved to {}'.format(outfn))

