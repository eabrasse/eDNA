#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
solve for population density and displacement kernel
"""

# setup
import numpy as np
from scipy import optimize
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import cmocean as cmo
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import random
from matplotlib.collections import LineCollection

tab20c = plt.get_cmap('tab20c',20)
plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'

#define spatial grid
x = np.linspace(-10,10,500)

#define some coefficients
sigma =0.5
mean = 3.5

#laplace
# pdf = 0.5*c0 * np.exp(-c0*np.abs(x))
#gaussian
pdf = (1/(sigma*np.sqrt(2*np.pi))) * np.exp(-0.5*((x-mean)/sigma)**2)

#normalize
pdf = pdf/np.sum(pdf)

#calculate cumulative density function
cdf = np.cumsum(pdf)

# samples
nsamples = 100
ntests = 25
samples = np.zeros((ntests,nsamples))
running_mean = np.zeros((ntests,nsamples))
running_var = np.zeros((ntests,nsamples))
randnum = np.random.rand(ntests,nsamples)
for i in range(nsamples):
    for t in range(ntests):
        ind = np.argmin(np.abs(randnum[t,i]-cdf))
        samples[t,i] = x[ind]
    running_mean[:,i] = np.mean(samples[:,:i],axis=1)
    running_var[:,i] = np.var(samples[:,:i],axis=1)



fig,axs = plt.subplots(figsize=(12,8),nrows=2,ncols=3)
#plot pdf
axs[0,0].plot(x,pdf,color='goldenrod',linestyle='solid')
axs[0,0].set_xlabel('x')
axs[0,0].set_ylabel('PDF')
axs[0,0].set_title('Gaussian Distribution PDF')

#plot cdf
axs[0,1].plot(x,cdf,color='green',linestyle='solid')
axs[0,1].set_xlabel('x')
axs[0,1].set_ylabel('CDF')
axs[0,1].set_title('Gaussian Distribution CDF')

#plot histogram of samples
counts, bins = np.histogram(samples[0,:],bins=np.arange(x.min(),x.max(),0.2))
axs[1,0].scatter(bins[1:],counts,facecolor='cornflowerblue',edgecolor='k')
axs[1,0].set_xlabel('x')
axs[1,0].set_ylabel('histogram of samples')
axs[1,0].set_title('Example sample distribution')

#plot running mean & variance
lw = 0.5
alp = 0.2
axs[1,1].axhline(y=mean,lw=1.,linestyle='dashed',color=tab20c(0),label='real mean')
axs[1,1].plot(range(nsamples),running_mean[0,:],label='running mean',alpha=alp,lw=lw,color=tab20c(0))
axs[1,1].axhline(y=sigma**2,lw=1.,linestyle='dashed',color=tab20c(4),label='real var')
axs[1,1].plot(range(nsamples),running_var[0,:],label='running var',alpha=alp,lw=lw,color=tab20c(4))
for t in range(1,ntests):
    axs[1,1].plot(range(nsamples),running_mean[t,:],label=None,alpha=alp,lw=lw,color=tab20c(0))
    axs[1,1].plot(range(nsamples),running_var[t,:],label=None,alpha=alp,lw=lw,color=tab20c(4))
axs[1,1].set_xlabel('Number of samples')
axs[1,1].legend()
axs[1,1].set_title('Mean and variance as function of samples')

# plot change in running mean
# running_mean_diff = running_mean[:,1:]-running_mean[:,:-1]
segs = np.zeros((ntests,nsamples-1,2))
running_mean_stag = 0.5*(running_mean[:,1:]+running_mean[:,:-1])
perc_change_running_mean = (running_mean[:,1:]-running_mean[:,:-1])/running_mean_stag
segs[:,:,1]=perc_change_running_mean
segs[:,:,0]=np.tile(np.reshape(np.arange(1,nsamples),(1,nsamples-1)),(ntests,1))
line_segments = LineCollection(segs,alpha=alp,lw=lw,color=tab20c(0))
axs[0,2].add_collection(line_segments)
axs[0,2].set_xlabel('Number of samples')
axs[0,2].set_ylabel('Percent change in running mean')
axs[0,2].autoscale()
axs[0,2].axhline(0.01,linestyle='dashed',color=tab20c(0))
axs[0,2].axhline(-0.01,linestyle='dashed',color=tab20c(0))


first_sample_below_1perc = np.zeros((ntests))
for t in range(ntests):
    first_sample_below_1perc[t] = np.argwhere(np.abs(perc_change_running_mean[t,:])<0.01)[0]
axs[1,2].scatter(range(ntests),first_sample_below_1perc,c=tab20c(1))
axs[1,2].set_xlabel('test')
axs[1,2].set_ylabel('first sample below 1 percent')
axs[1,2].axhline(np.mean(first_sample_below_1perc),linestyle='dashed',color=tab20c(0))


plt.subplots_adjust(left=0.07,right=0.9,bottom=0.1,top=0.95,hspace=0.3,wspace=0.2)
plt.show(block=False)
plt.pause(0.1)
# outfn = home + 'plots/displacement_kernel/figure_{:03d}.png'.format(tt)
# plt.savefig(outfn)
# plt.close()
