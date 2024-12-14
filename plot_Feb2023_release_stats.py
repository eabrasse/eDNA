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
paired = plt.get_cmap('Paired',12)
tz = pytz.utc
pdt = pytz.timezone('America/Vancouver')

# load in data
home = '/Users/elizabethbrasseale/Projects/eDNA/'


data_fn = home+'data/Feb2023_hist_counts.p'
D = pickle.load(open(data_fn,'rb'))
for key in D.keys():
    locals()[key] = D[key]

fw,fh = efun.gen_plot_props()
fig = plt.figure(figsize=(3*fw,2*fh))

#first make sure totals add up
ax = fig.add_subplot(2,1,1)

hist_sum = np.sum(np.sum(hist_particles_bin,2),1)
ax.plot(dt_list0,total_released_particles,color='k',label='tot rel')
ax.plot(dt_list0,hist_sum,color=paired(1),label=r'$\Sigma$ hist')
ax.plot(dt_list0,surface_particles,color=paired(0),linestyle='dashed',label='surf')
ax.plot(dt_list0,outfield_particles,color=paired(3),label='outfield')
ax.plot(dt_list0,deep_particles,color=paired(5),label='deep',ls='dashed')
ax.plot(dt_list0,hist_sum+deep_particles+outfield_particles,color='gray',linestyle='dashed',label=r'$\Sigma$ hist + outfield + deep')
ax.set_ylabel('# particles')
ax.set_xlabel('Date time')
ax.legend()

# now look at change across time steps
dt = dt_list0[1]-dt_list0[0]
dt_list1 = [dt_list0[i]+dt for i in range(len(dt_list0))[:-1]]

Dprod = np.diff(total_released_particles)
Dpart = np.diff(hist_sum)
Dsink = np.diff(deep_particles)
Dleave = np.diff(outfield_particles)

ax1 = fig.add_subplot(2,1,2)
ax1.plot(dt_list1[:-1],Dprod[:-1],color='k',label='prod')
ax1.plot(dt_list1[:-1],Dpart[:-1],color=paired(1),label=r'$\frac{D}{Dt}\Sigma$ hist')
ax1.plot(dt_list1[:-1],Dsink[:-1],color=paired(5),linestyle='dashed',label='sink')
ax1.plot(dt_list1[:-1],Dleave[:-1],color=paired(3),label='left domain')
ax1.plot(dt_list1[:-1],Dprod[:-1]-Dsink[:-1]-Dleave[:-1],color=paired(0),linestyle='dashed',label='Prod-Sink-Left')

ax1.set_ylabel(r'$\Delta \#$ particles per hour')
ax1.set_xlabel('Date time')
ax1.legend()



# fig.subplots_adjust(left=0.2,right=0.98,bottom=0.2)
plt.show(block=False)
plt.pause(0.1)
