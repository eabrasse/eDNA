"""
Plot results of a particle tracking experiment.
"""


import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import matplotlib

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()

paired = plt.get_cmap('Paired',12)
myblue = [p for p in paired(5)]
colors = [(myblue[0], myblue[1], myblue[2],c) for c in np.linspace(0,1,100)]
newBlues = matplotlib.colors.LinearSegmentedColormap.from_list('mycmap', colors, N=100)

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
hh = dsg.h.values
maskr = dsg.mask_rho.values
#

# subsample output for plotting
npmax = 300 # max number of points to plot
step = max(1,int(np.floor(NP/npmax)))

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

# PLOTTING
plt.close('all')
pfun.start_plot(figsize=(14,8))
fig = plt.figure()

# MAP
# set domain limits
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
    
# ax = fig.add_subplot(121)
ax = fig.gca()
zm = -np.ma.masked_where(maskr==0, hh)
# plt.pcolormesh(lonp, latp, zm, vmin=-100, vmax=0,
    # cmap='terrain', alpha=.25)
pfun.add_coast(ax)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title(exp_name.strip('/'))
# add the tracks (packed [time, particle])
# regular spaghetti plots
# ax.plot(lon, lat, '-k', linewidth=.2)
ax.plot(lon[0,0], lat[0,0], marker='*',mec='k',mfc='yellow', markersize=10,alpha=1,zorder=500)
# ax.plot(lon[-1,:], lat[-1,:], 'or', alpha=.3)
nbins = 50
bin_lon_edges=np.linspace(aa[0], aa[1],nbins+1)
bin_lat_edges=np.linspace(aa[2], aa[3],nbins+1)
xx, yy = np.meshgrid(bin_lon_edges[:-1]+0.5*(bin_lon_edges[1]-bin_lon_edges[0]),bin_lat_edges[:-1]+0.5*(bin_lat_edges[1]-bin_lat_edges[0]))
hist = np.histogram2d(lon[-1,:],lat[-1,:],bins=[bin_lon_edges,bin_lat_edges])
vmax=np.percentile(hist[0][hist[0]>0],90)
p=ax.pcolormesh(xx,yy,hist[0].transpose(),vmax=vmax,cmap=newBlues)
plt.colorbar(p)

# # time series
# td = (ot_vec - ot_vec[0])/86400
# tv_list = ['z', 'cs', 'salt', 'temp']
# #tv_list = ['u', 'v', 'lon', 'lat']
# ntv = len(tv_list)
# for ii in range(ntv):
#     tv = tv_list[ii]
#     NC = 2
#     ax = fig.add_subplot(ntv,NC, (ii+1)*NC)
#     v = dsr[tv].values[:,::step]
#     v[~ib_mask] = np.nan
#     ax.plot(td, v, lw=.5, alpha=.2)
#     ax.text(.05, .05, tv, fontweight='bold', transform=ax.transAxes)
#     if ii == ntv-1:
#         ax.set_xlabel('Time (days)')

#plt.show()
outfn = '/data2/pmr4/eab32/etools/plots/hc_dolph_surf_sh1_10k_heatmap_tight.png'
plt.savefig(outfn)
pfun.end_plot()

print('saved to {}'.format(outfn))
dsr.close()
dsg.close()

