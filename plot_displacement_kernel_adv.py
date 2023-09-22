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
from scipy import special


plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'
tab20c = plt.get_cmap('tab20c',20)
vmin  = 0
vmax = 0.5

#define spatial grid
x = np.linspace(-10,10,101)

#define some coefficients
D = 1
alpha = 1
sigma = np.sqrt(2*D/alpha)
b = np.sqrt(D/alpha)
v=1

# solve for steady state analytical displacement kernel
c0 = 1/b
c1 = 1/(np.pi*(sigma**2))
c2 = 2/(sigma**2)

n= np.arange(0,10)
k1D = np.zeros((len(n),len(x)))

a1 = (v+np.sqrt(v**2+4*alpha*D))/(2*D)
a2 = (v-np.sqrt(v**2+4*alpha*D))/(2*D)

for nn in n:
    k1D[nn,x>=0] = 0.5*(nn+1)*c0 * np.exp(-a1*x[x>=0])
    k1D[nn,x<0] = 0.5*(nn+1)*c0 * np.exp(-a2*x[x<0])


samples = 1
fig,axs = plt.subplots(nrows=1,ncols=2,figsize=(11,5))
rainbow = plt.get_cmap('rainbow',len(n))

n_real = 4
x_start = -4
x_step = 1
k_uncertainty = 0.1

x_pos = [] #position where sampling
k_samp = [] #value of sample

x_pos.append(x_start)
ax = axs[0]
ax1 = axs[1]
for s in range(samples):
    #find index of x position
    xi = np.argmin(np.abs(x-x_pos[-1]))
    # 'take' sample
    k_samp.append(k1D[n_real,xi])
    ax.axhline(k_samp[-1],linestyle='dashed',color='k')
    # loop through possible number of animals
    for nn in n:
        # find where the sample is within uncertainty - it's likely a segment of k
        ki = np.argwhere(np.abs(k_samp[-1]-k1D[nn,:])<k_uncertainty)
        ki0 = 1+np.argwhere((ki[1:]-ki[:-1])>1)
        ki0 = [ki[t+1] for t in range(len(ki)-1) if (ki[t+1]-ki[t])>1]
        ki1 = [ki[t] for t in range(len(ki)-1) if (ki[t+1]-ki[t])>1]
        ki0 = np.insert(ki0,0,ki[0])
        ki1 = np.append(ki1,ki[-1])
        
        # count number of segments
        if len(ki0)==len(ki1):
            klen = len(ki0)
        else:
            klen = int(np.min([len(ki0),len(ki1)]))
        # loop through segments where k is within uncertainty of k_samp
        for kseg in range(klen):
            # read in axes that are beginning and end of segment
            k0 = ki0[kseg]
            k1 = ki1[kseg]
            # find most likely position given segment
            kii = k0+np.argmin(np.abs(k_samp[-1]-k1D[nn,k0:k1]))
            
            xerr = np.reshape(np.array([x[kii]-x[k0],x[k1]-x[kii]]),(2,1))
            # ax.plot(x-x[kii],k1D[nn,:],color=rainbow(nn),linestyle='solid')
            ax.fill_between(x-x[kii],k1D[nn,:]-k_uncertainty,k1D[nn,:]+k_uncertainty,color=rainbow(nn),alpha=0.3)
            ax1.errorbar(-x[kii],nn+1,xerr=xerr,linestyle='solid',color=rainbow(nn),marker='o',mfc=rainbow(nn),mec='k')

# ax.axhline(k_sam,color='k',linestyle='dashed')
ylim = [-0.3,5.3]
xlim = [-6,6]
ax.set_ylim(ylim)
ax.set_xlim(xlim)
ax.set_xlabel('Distance from you (m)')
ax.set_ylabel('eDNA conc')
# ax.set_title('Location and number of dolphins')

xtext = 0.1
ytext = 0.9
ax.text(xtext,ytext,'Number of dolphins',color='k',transform=ax.transAxes)
epsilon = 0.04
for nn in n:
    ax.text(xtext,ytext-(nn+1)*epsilon,'{}'.format(nn+1),color=rainbow(nn),transform = ax.transAxes)

ax.axvline(0,linestyle='dashed',color='k')
ax.text(-0.03,np.mean(ylim),'your location',ha='right',va='center',rotation=90)
ax.text(xlim[0]+0.15*(xlim[1]-xlim[0]),k_samp[-1]+0.03,'eDNA measured')
ax.text(0.7,0.65,'viable distribution functions\nwhich reproduce measured eDNA',rotation = 70,va='bottom')

ax1.set_title('Possible number and locations of dolphins producing eDNA')
ax1.set_xlabel('Distance from you (m)')
ax1.set_ylabel('Number of dolphins')

plt.subplots_adjust(left=0.07,right=0.9,bottom=0.1,top=0.95,hspace=0.3,wspace=0.2)
plt.show(block=False)
plt.pause(0.1)
