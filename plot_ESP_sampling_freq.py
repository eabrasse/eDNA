#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
solve for population density and displacement kernel
"""

# setup
import numpy as np
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec
import itertools
import seaborn as sns
import matplotlib.dates as mdates

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'
tab20c = plt.get_cmap('tab20c',20)

# load in samples data
ESP_fn = home+'data/03_ESP1_Feb2023_hourly.csv'
df = pd.read_csv(ESP_fn,sep=',',engine='python')
# df['datetime']=pd.to_datetime(df.date)

DNA_data = df['PB_quantity_mean'].dropna().values

fig, axs = plt.subplots(1,4,figsize=(16,5))
df['datetime']=pd.to_datetime(df.date)
axs[0].scatter(df.datetime,df.PB_quantity_mean,color=tab20c(0))#edgecolors=tab10(0),c='None')
axs[0].set_yscale('log')
axs[0].grid()
axs[0].set_ylabel(r'Conc (copies/$\mathrm{\mu L}$)')
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%b %-d %-H:%m"))
plt.setp( axs[0].xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
axs[0].plot(df.datetime[df.PB_quantity_mean==df.PB_quantity_mean.max()],df.PB_quantity_mean[df.PB_quantity_mean==df.PB_quantity_mean.max()],marker='*',mec='k',mfc=tab20c(0),markersize=15)
axs[0].annotate('max',xy=(df.datetime[df.PB_quantity_mean==df.PB_quantity_mean.max()],df.PB_quantity_mean[df.PB_quantity_mean==df.PB_quantity_mean.max()]),xytext=(-10,0.),
    textcoords='offset points',color=tab20c(0),va='center',ha='right')
ylim = axs[0].get_ylim()
axs[0].set_ylim((ylim[0],ylim[1]*(1+np.sqrt(2))))
nt = len(DNA_data)
D = {}
for Nsample in range(1,nt+1):
    D[Nsample] = {}
    D[Nsample]['mean'] = []
    D[Nsample]['mean_wo_max'] = []
DNA_data_max = DNA_data.max()
# [[[k for k in range(i,nt,step+1)] for step in range(nt-i)] for i in range(nt)] prints all desired combos
combs = [[[range(start,stop,step) for step in range(1,1+stop-start)] for stop in range(start+1,nt)] for start in range(nt)]
for start_combs in combs: # will be a list of ranges of inds
    for startstop_combs in start_combs:
        for startstopstep_combs in startstop_combs:
            Nsample = len(startstopstep_combs)
            if Nsample>0:
                mean = np.mean(DNA_data[startstopstep_combs])
                D[Nsample]['mean'].extend([mean])
                if DNA_data_max not in DNA_data[startstopstep_combs]:
                    D[Nsample]['mean_wo_max'].extend([mean])

stddev = np.zeros((nt-1))
rmse = np.zeros((nt-1))
nrmse = np.zeros((nt-1))
stddev_wo_max = np.zeros((nt-1))
rmse_wo_max = np.zeros((nt-1))
nrmse_wo_max = np.zeros((nt-1))

Nsamples = np.array(range(1,nt+1))
keys_wo_max = [Nsample for Nsample in Nsamples if len(D[Nsample]['mean_wo_max'])>0]
maxkey_wo_max = max(keys_wo_max)

for Nsample in Nsamples:
    y = D[Nsample]['mean']
    ny = len(y)
    x = [Nsample]*ny
    
    y_wo_max = D[Nsample]['mean_wo_max']
    ny_wo_max = len(y_wo_max)
    x_wo_max = [Nsample]*ny_wo_max
    
    y_w_max = [y0 for y0 in y if y0 not in y_wo_max]
    ny_w_max = len(y_w_max)
    x_w_max = [Nsample]*ny_w_max
    axs[1].scatter(x_wo_max,y_wo_max,color=tab20c(0),alpha=0.15,s=[10],zorder=20)
    axs[1].scatter(x_w_max,y_w_max,color=tab20c(0),alpha=0.15,s=[10],zorder=20)
        
    if Nsample<Nsamples.max(): #cut off zero to plot in log space
        stddev[Nsample-1] = np.std(y)
        rmse[Nsample-1] = np.sqrt((np.sum((y-D[41]['mean'][0])**2))/ny)
        nrmse[Nsample-1] = np.sqrt((np.sum((y-D[41]['mean'][0])**2))/D[41]['mean'][0]**2)
        

        
        stddev_wo_max[Nsample-1] = np.std(y_wo_max)
        rmse_wo_max[Nsample-1] = np.sqrt((np.sum((y_wo_max-D[maxkey_wo_max]['mean_wo_max'][0])**2))/ny_wo_max)
        nrmse_wo_max[Nsample-1] = np.sqrt((np.sum((y_wo_max-D[maxkey_wo_max]['mean_wo_max'][0])**2))/D[maxkey_wo_max]['mean_wo_max'][0]**2)

# axs[0].set_xlabel('Number of samples')
axs[1].set_ylabel('Mean')
axs[1].text(17,11428,'subsets containing max',color=tab20c(0),ha='center',va='center',rotation=-20)
axs[1].text(22.5,500,'subsets not\ncontaining max',color=tab20c(0),ha='center',va='center',rotation=0)


axs[2].plot(Nsamples[:-1],stddev,color=tab20c(4))
axs[2].plot(Nsamples[:-1],stddev_wo_max,color=tab20c(4),linestyle='dashed')
# axs[1].set_xlabel('Number of samples')
axs[2].set_ylabel('Standard deviation')
axs[2].text(13,3880,'all subsets',color=tab20c(4),ha='center',va='center',rotation=-20)
axs[2].text(7,480,'only subsets\nw/o max',color=tab20c(4),ha='center',va='center',rotation=-30)


axs[3].plot(Nsamples[:-1],nrmse,color=tab20c(8))
axs[3].plot(Nsamples[:-1],nrmse_wo_max,color=tab20c(8),linestyle='dashed')
# axs[2].set_xlabel('Number of samples')
axs[3].set_ylabel('NRMSE')

for ax in axs[1:]:
    ax.set_xlabel('Number of samples')
    ax.set_yscale('log')
    ax.grid()
# next violin plot
# sns.violinplot( x=Nsample, y=mean_list,ax=axs[1],color=tab20c(4))
# axs[1].violinplot( data = mean_list, pos=Nsample,color=tab20c(4))

axs[0].text(0.1,0.96,'a) ESP samples near net pen',transform=axs[0].transAxes,ha='left',va='top')
axs[1].text(0.1,0.96,'b) Means using subset of samples\n(vary sampling start, duration, & freq)',transform=axs[1].transAxes,ha='left',va='top')
axs[2].text(0.1,0.96,'c) Standard deviation of possible means',transform=axs[2].transAxes,ha='left',va='top')
axs[3].text(0.1,0.96,'d) Normalized root mean square error\n of possible means',transform=axs[3].transAxes,ha='left',va='top')

plt.subplots_adjust(left=0.05,right=0.99,bottom=0.15,top=0.95,wspace=0.2)
plt.show(block=False)
plt.pause(0.1)