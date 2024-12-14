import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import pandas as pd
import string
import numpy as np
import pytz
import pickle
from matplotlib.gridspec import GridSpec
import calendar
import efun
from scipy.stats import linregress
import matplotlib.dates as mdates
from matplotlib.dates import DayLocator, HourLocator, MonthLocator,DateFormatter, drange

plt.close('all')
tab10 = plt.get_cmap('tab10',10)
tz = pytz.utc
pdt = pytz.timezone('America/Vancouver')

# load in data
home = '/Users/elizabethbrasseale/Projects/eDNA/'


data_fn = home+'data/Feb2023_surfaceparticlepositionstats.p'
D = pickle.load(open(data_fn,'rb'))


fw,fh = efun.gen_plot_props()
fig = plt.figure(figsize=(2*fw,2*fh))
ax_avg_dist = fig.add_subplot(3,1,1)
ax_dens = fig.add_subplot(3,1,2)
ax_z = fig.add_subplot(3,1,3)
tscale = 1/3600
all_release_depths = {}
all_release_densities = {}
all_release_distances = {}
t_list = np.arange(len(D.keys()))
for t in t_list:
    all_release_depths[t] = []
    all_release_densities[t] = []
    all_release_distances[t] = []
count=0
for key in D.keys():
    rel = D[key]
    
    density = 1/(2*rel['avg_dist_from_COM']**2)
    ax_avg_dist.plot(tscale*rel['t'],rel['avg_dist_from_COM'],color=tab10(0),alpha=0.1)
    ax_dens.plot(tscale*rel['t'],density,color=tab10(2),alpha=0.1)
    ax_z.plot(tscale*rel['t'],rel['deep_masked_count'],color=tab10(1),alpha=0.1)
    
    for t in tscale*rel['t']:
        all_release_depths[t].append(rel['deep_masked_count'][int(t)])
        all_release_densities[t].append(density[int(t)])
        all_release_distances[t].append(rel['avg_dist_from_COM'][int(t)])



# if I expect them to decay exponentially...
tco = 24 # tco = time cut off
mean_release_density = [np.mean(all_release_densities[t]) for t in t_list]
slope, intercept, r_value, p_value, std_err = linregress(t_list[:tco],np.log(mean_release_density[:tco]))
yy = np.exp(intercept)*np.exp(t_list[:tco]*slope)
ax_dens.plot(t_list,mean_release_density,color='k',linestyle='dashed',lw=2)
ax_dens.plot(t_list[:tco],yy,color='m',linestyle='dashed',lw=2)
ax_dens.text(0.5,0.5,r'$k_d=$'+f'{slope:.3} '+r'$\mathrm{hr}^{-1}$',transform=ax_dens.transAxes,color='m')

mean_release_depth = [np.mean(all_release_depths[t])+100 for t in t_list]
slope, intercept, r_value, p_value, std_err = linregress(t_list[:],np.log(mean_release_depth[:]))
yy = np.exp(intercept)*np.exp(t_list[:]*slope)
ax_z.plot(t_list,np.array(mean_release_depth)-100,color='k',linestyle='dashed',lw=2)
ax_z.plot(t_list[:],np.array(yy)-100,color='m',linestyle='dashed',lw=2)
ax_z.text(0.1,0.1,r'$k_s=$'+f'{slope:.3} '+r'$\mathrm{hr}^{-1}$',transform=ax_z.transAxes,color='m')


mean_release_distance = [np.mean(all_release_distances[t])+100 for t in t_list]
ax_avg_dist.plot(t_list,mean_release_distance,color='k',linestyle='dashed',lw=2)
ax_avg_dist.set_ylabel('Distance from COM')
ax_avg_dist.set_xlabel('')
# ax_avg_dist.set_yscale('log')
ax_avg_dist.text(0.1,0.9,'Average distance of surface particles to COM',transform=ax_avg_dist.transAxes)


ax_dens.set_ylabel('1/Distance from COM cubed (m-2)')
ax_dens.set_xlabel('')
ax_dens.set_yscale('log')
ax_dens.text(0.1,0.9,'Average surface particle concentration',transform=ax_dens.transAxes)

# > ax_dist.set_yscale('log')
ax_z.set_ylabel('N particles')
ax_z.set_xlabel('Time since release (hours)')
ax_z.text(0.1,0.9,'Total depth masked particles',transform=ax_z.transAxes)

# for ax in ax_avg_dist,ax_z:
#     ax.set_xlim([0,24])

fig.subplots_adjust(top=0.98,left=0.2,right=0.98,bottom=0.1)
plt.show(block=False)
plt.pause(0.1)
