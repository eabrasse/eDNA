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

atoz = string.ascii_lowercase
Ldir = Lfun.Lstart()

tab10 = plt.get_cmap('tab10',10)

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

nbins = 100
bin_lon_edges=np.linspace(aa[0], aa[1],nbins+1)
bin_lat_edges=np.linspace(aa[2], aa[3],nbins+1)
xx, yy = np.meshgrid(bin_lon_edges[:-1]+0.5*(bin_lon_edges[1]-bin_lon_edges[0]),bin_lat_edges[:-1]+0.5*(bin_lat_edges[1]-bin_lat_edges[0]))


tvary_denom = 0
pulse_hist_decay = np.zeros((len(release_list),np.shape(xx)[0],np.shape(xx)[1]))
tvary_hist_decay = np.zeros(np.shape(xx))
deltaT_list = []
count=0
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
    if release==release_list[0]:
        NP = lon.shape[0]

    pulse_hist_decay[count,:] = rel['decay']*hist[0].T/NP
    tvary_hist_decay += rel['decay']*rel['C0']*hist[0].T
    
    tvary_denom += rel['C0']*NP
    
    deltaT_list.append(rel['deltaT']/3600)
    
    count+=1
dT_inds = np.argsort(deltaT_list[::-1])
pulse_hist_decay=pulse_hist_decay[dT_inds,:]

#normalize
tvary_hist_decay=tvary_hist_decay/tvary_denom

# sort indexes
tv_inds = tvary_hist_decay.flatten().argsort()
xsort = 100*tvary_hist_decay.flatten()[tv_inds]
xmean = tvary_hist_decay.mean()

nrel = len(release_list)
rmse = {'values':np.zeros((nrel)),'label':'RMSE','best':0}
nrmse = {'values':np.zeros((nrel)),'label':'NRMSE','best':0}
slope = {'values':np.zeros((nrel)),'label':'Slope of\nbest fit line','best':1}
intercept = {'values':np.zeros((nrel)),'label':'Intercept of best fit line','best':0}
wss = {'values':np.zeros((nrel)),'label':'WSS','best':1}
Rsquared = {'values':np.zeros((nrel)),'label':r'$\mathrm{R}^{2}$','best':1}

for i in range(nrel):
    rmse['values'][i] = np.sqrt((np.sum((tvary_hist_decay.flatten()-pulse_hist_decay[i,:].flatten())**2))/len(tvary_hist_decay.flatten()))
    nrmse['values'][i] = rmse['values'][i]/xmean

    ysort = [100*pulse_hist_decay[i,:].flatten()[tv_ind] for tv_ind in tv_inds]
    slope['values'][i],intercept['values'][i] = np.polyfit(xsort,ysort,1)
    # xfit = np.linspace(0,xymax,10)
    # yfit = pfit[0]*xfit +pfit[1]
    # ax2.plot(xfit,yfit,linestyle='dashed',color='magenta',zorder=26)
    R = np.corrcoef(xsort,ysort)[0,1]
    Rsquared['values'][i] = R**2
    wss['values'][i] = efun.willmott(xsort,ysort)

metrics2plot = [nrmse,slope,wss,Rsquared]

# PLOTTING
plt.close('all')
fw, fh = efun.gen_plot_props()
figsize = (fw*2.5,fh)
fig,axs = plt.subplots(1,len(metrics2plot),figsize=figsize)
deltaT_list.sort(reverse=True)
for met in range(len(metrics2plot)):
    metric = metrics2plot[met]
    ax = axs[met]
    
    ax.plot(deltaT_list,metric['values'],color=tab10(met),zorder=15)
    ax.text(0.1,0.9,'{}) {}'.format(atoz[met],metric['label']),color='k',transform=ax.transAxes,zorder=100)
    
    best_ind = np.argmin(np.abs(metric['values']-metric['best']))
    ax.plot(deltaT_list[best_ind],metric['values'][best_ind],marker='*',linestyle='none',markersize=10,mec='k',mfc=tab10(met),zorder=50)
    ax.axvline(deltaT_list[best_ind],color='k',linestyle='dashed',zorder=12)
    ax.text(0.5,0.5,'best fit release =\n{} hours before snapshot'.format(deltaT_list[best_ind]),transform=ax.transAxes,color='k',zorder=70,ha='center',fontsize=8)
    
    ax.set_xlabel('Release #')


#plt.show()
fig.subplots_adjust(bottom=0.08,top=0.98,left=0.08,right=0.98,wspace=0.2)
outfn = '/data2/pmr4/eab32/etools/plots/hc_dolph_3d_compare_pulse_releases.png'
fig.savefig(outfn)
print('saved to {}'.format(outfn))


plt.close('all')
