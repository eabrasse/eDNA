"""
Plot results of a particle tracking experiment.
"""


import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import matplotlib
import matplotlib.dates as mdates
import pickle
from lo_tools import Lfun,zfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()
from datetime import timedelta, datetime

paired = plt.get_cmap('Paired',12)
myblue = [p for p in paired(5)]
colors = [(myblue[0], myblue[1], myblue[2],c) for c in np.linspace(0,1,100)]
newBlues = matplotlib.colors.LinearSegmentedColormap.from_list('mycmap', colors, N=100)

# useful plotting tools, colormap and test box props
landcol = 'lightgray'
seacol = 'white'

cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])


# identify directory for output
outdir = '/data2/pmr4/eab32/etools/plots/hc_dolph_particles_tight_mov/'

# mask path
mask_fn = '/data2/pmr4/eab32/LO_output/extract/hc11_v1_uu0k/wetdry_mask_tides.p'
Dmask = pickle.load(open(mask_fn,'rb'))
wd_mask = Dmask['wetdry_mask_rho'][:]
wd_dt_list = Dmask['dt_list'][:]

# Choose an experiment and release to plot.
in_dir0 = Ldir['LOo'] / 'tracks'
exp_name = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
    itext='** Choose experiment from list **', last=False)
rel = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)

# get Datasets
fn = in_dir0 / exp_name / rel
fng = in_dir0 / exp_name / 'grid.nc'
dsr = xr.open_dataset(fn, decode_times=False)
dsg = xr.open_dataset(fng)

NT, NP = dsr.lon.shape

# get a list of datetimes
ot_vec = dsr.ot.values
dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]

# gather some fields, for convenience
lonp, latp = pfun.get_plon_plat(dsg.lon_rho.values, dsg.lat_rho.values)
lonr = dsg.lon_rho.values
latr = dsg.lat_rho.values
hh = dsg.h.values
maskr = dsg.mask_rho.values

# get tidal signal
lonvec = lonr[0,:]
latvec = latr[:,0]
ks0201_lon = -122.73
ks0201_lat = 47.75833
i_ks0201 = zfun.find_nearest_ind(lonvec,ks0201_lon)
j_ks0201 = zfun.find_nearest_ind(latvec,ks0201_lat)
zeta = Dmask['zeta'][:,j_ks0201,i_ks0201]



# subsample output for plotting
# npmax = 300 # max number of points to plot
# step = max(1,int(np.floor(NP/npmax)))

# lon = dsr.lon.values[:,::step]
# lat = dsr.lat.values[:,::step]
lon = dsr.lon.values[:]
lat = dsr.lat.values[:]

# make a mask that is False from the time a particle first leaves the domain
# and onwards
AA = [dsg.lon_rho.values[0,0], dsg.lon_rho.values[0,-1],
        dsg.lat_rho.values[0,0], dsg.lat_rho.values[-1,0]]
ib_mask = np.ones(lon.shape, dtype=bool)
ib_mask[lon < AA[0]] = False
ib_mask[lon > AA[1]] = False
ib_mask[lat < AA[2]] = False
ib_mask[lat > AA[3]] = False
NTS, NPS = lon.shape
for pp in range(NPS):
    tt = np.argwhere(ib_mask[:,pp]==False)
    if len(tt) > 0:
        ib_mask[tt[0][0]:, pp] = False

# and apply the mask to lon and lat
lon[~ib_mask] = np.nan
lat[~ib_mask] = np.nan

#generate spatial histogram bins 
# define domain limits
if False:
    # plot full domain
    aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
else:
    # automatically plot region of particles, with padding
    pad = .01
    # aa = [np.nanmin(lon) - pad, np.nanmax(lon) + pad,
    # np.nanmin(lat) - pad, np.nanmax(lat) + pad]
    
    # plot region around delta pier, using release point as indicator
    aa = [lon[0,0]-pad, lon[0,0]+pad,
       lat[0,0]-pad, lat[0,0]+pad]
nbins = 50
bin_lon_edges=np.linspace(aa[0], aa[1],nbins+1)
bin_lat_edges=np.linspace(aa[2], aa[3],nbins+1)
xx, yy = np.meshgrid(bin_lon_edges[:-1]+0.5*(bin_lon_edges[1]-bin_lon_edges[0]),bin_lat_edges[:-1]+0.5*(bin_lat_edges[1]-bin_lat_edges[0]))
vmax=20 # count cap for colormap

# get number of time steps to loop over
nt, ni = lon.shape

for tt in range(nt):
    # PLOTTING
    plt.close('all')
    pfun.start_plot(figsize=(10,8))
    fig = plt.figure()
    
    ax = fig.add_subplot(121)
    # ax = fig.gca()
    zm = -np.ma.masked_where(maskr==0, hh)
    # plt.pcolormesh(lonp, latp, zm, vmin=-100, vmax=0,
        # cmap='terrain', alpha=.25)
        
    #add mask and wetdry mask
    wd_td_list = [np.abs(wd_dt - dt_list[tt]) for wd_dt in wd_dt_list]
    wd_ind = wd_td_list.index(min(wd_td_list))
    ax.pcolormesh(lonp,latp,maskr,cmap=cmap_mask,shading='flat',zorder=0)
    ax.contour(lonr,latr,wd_mask[wd_ind,:,:],levels=[0.5],colors=['k'],linewidths=[1.5])
    
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(exp_name.strip('/'))
    # add the tracks (packed [time, particle])
    # regular spaghetti plots
    # ax.plot(lon, lat, '-k', linewidth=.2)
    ax.plot(lon[0,0], lat[0,0], marker='*',mec='k',mfc='yellow', markersize=10,alpha=1,zorder=500)
    # ax.plot(ks0201_lon, ks0201_lat, marker='d',mec='k',mfc='orange', markersize=10,alpha=1,zorder=500) # note: outside of axes
    # ax.plot(lon[-1,:], lat[-1,:], 'or', alpha=.3)

    # hist = np.histogram2d(lon[tt,:],lat[tt,:],bins=[bin_lon_edges,bin_lat_edges])
    # # vmax=np.percentile(hist[0][hist[0]>0],90)
    # p=ax.pcolormesh(xx,yy,hist[0].transpose(),vmax=vmax,cmap=newBlues)
    # plt.colorbar(p)
    ax.scatter(lon[tt,:],lat[tt,:],c='green',alpha=0.3)

    # time series
    ax = fig.add_subplot(122)
    ax.plot(wd_dt_list,zeta,color='gray')
    ax.plot(wd_dt_list[wd_ind],zeta[wd_ind],marker='o',mfc=paired(1))
    ax.set_xlabel('2023 Date, time (UTC)')
    ax.set_xlim([datetime(2023,1,30,22),datetime(2023,2,2,2)])
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d %H:%M"))
    plt.setp( ax.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
    ax.set_ylabel('Tide height (m)')
    ax.set_aspect(0.4)
    ax.text(0.1,0.9,'ks0201',transform=ax.transAxes)

    outfn = outdir+f'fig_{tt:03}.png'
    plt.savefig(outfn)
    pfun.end_plot()

    print('saved to {}'.format(outfn))
dsr.close()
dsg.close()

