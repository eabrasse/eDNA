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
vmax = 1.0
ylim = [0,0.55]

#define spatial grid
x = np.linspace(-10,10,101)

#define some coefficients
D = 1
alpha = 0.5
sigma = np.sqrt(2*D/alpha)
b = np.sqrt(D/alpha)

# now solve for time dependent spin up
# initialize time grid
t = np.linspace(0,5,1000)

#solve for dx, dt
dx = x[1]-x[0]
dt = t[1]-t[0]

# initialize a reversing velocity
v = 2.0*np.sin(t*2*np.pi/(t[-1]-t[0]))
v0 = np.zeros(v.shape)

# initialize k solution
k_1noise_v = np.zeros((len(t),len(x)))
k_2noise_v = np.zeros((len(t),len(x)))

#production scale
prod = 10

#condition at t=0 (index 50 when x=np.linspace(-50,50,101))
production_1noise = np.zeros((len(t),len(x)))
production_1noise[:,50] = prod*(1+(np.cos(t*7*np.pi/(t[-1]-t[0]))+np.sin(t*2*np.pi/(t[-1]-t[0]))))
production_2noise = np.zeros((len(t),len(x)))
production_2noise[:,50] = 2*prod*(1+(np.cos(t*7*np.pi/(t[-1]-t[0]))+np.sin(t*2*np.pi/(t[-1]-t[0]))))
for production_noise in [production_1noise, production_2noise]:
    production_noise[production_noise<0]=0

k_1noise_v[0,:] += production_1noise[0,:]*dt
k_2noise_v[0,:] += production_2noise[0,:]*dt

for tt in range(1,len(t)):
    
    # noise, current
    for [k, production] in [[k_1noise_v,production_1noise],[k_2noise_v,production_2noise]]:
        dkdx_up = (k[tt-1,0:-2]-k[tt-1,1:-1])/dx
        dkdx_dn = (k[tt-1,1:-1]-k[tt-1,2:])/dx
        if v[tt-1]>=0:
            adv = v[tt-1]*dkdx_up
        elif v[tt-1]<0:
            adv = v[tt-1]*dkdx_dn
        diffusion = D*(dkdx_up-dkdx_dn)/dx
        decay = -k[tt-1,:]*alpha
        k[tt,1:-1] = k[tt-1,1:-1] + dt*(adv + diffusion + decay[1:-1]+ production[tt,1:-1])
    
    # # noise, current
    # dkdx_up_2noise_v = (k_2noise_v[tt-1,0:-2]-k_2noise_v[tt-1,1:-1])/dx
    # dkdx_dn_2noise_v = (k_2noise_v[tt-1,1:-1]-k_2noise_v[tt-1,2:])/dx
    # if v[tt-1]>=0:
    #     adv_2noise_v = v[tt-1]*dkdx_up_2noise_v
    # elif v[tt-1]<0:
    #     adv_2noise_v = v[tt-1]*dkdx_dn_n2oise_v
    # diffusion_2noise_v = D*(dkdx_up_2noise_v-dkdx_dn_2noise_v)/dx
    # decay_2noise_v = -k_2noise_v[tt-1,:]*alpha
    # k_2noise_v[tt,1:-1] = k_2noise_v[tt-1,1:-1] + dt*(adv_2noise_v + diffusion_2noise_v + decay_2noise_v[1:-1]+ production_2noise[tt,1:-1])
    

#plotting params
cmap = plt.get_cmap(cmo.cm.dense)
xx,tt = np.meshgrid(x,t)

#generate figure
fig,axs = plt.subplots(figsize=(10,8),nrows=2,ncols=2)

# plot all dye in time, space
for ax,k in [[axs[0,0],k_1noise_v],[axs[1,0],k_2noise_v]]:
    p=ax.pcolormesh(tt,xx,k,cmap=cmo.cm.matter,vmin=vmin,vmax=vmax)
    ax.axhline(0,linestyle='dashed',color='k')
    ax.set_aspect(0.25)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Distance (m)')

# add colorbar
cbaxes = inset_axes(axs[0,-1], width="4%", height="60%", loc='center right',bbox_transform=axs[0,-1].transAxes,bbox_to_anchor=(0.1,0.,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')

# include insets of production
pax0 = inset_axes(axs[0,0], width="15%", height="15%", loc='lower left',bbox_transform=axs[0,0].transAxes,bbox_to_anchor=(0.1,0.05,1,1))
pax1 = inset_axes(axs[1,0], width="15%", height="15%", loc='lower left',bbox_transform=axs[1,0].transAxes,bbox_to_anchor=(0.1,0.05,1,1))

pax0.plot(t,production_1noise[:,50],color=tab20c(0))
pax1.plot(t,production_2noise[:,50],color=tab20c(0))
    
for pax in [pax0,pax1]:
    pax.set_ylabel('eDNA production',fontsize=6)
    pax.set_ylim([-20,60])


# include insets with current
vax0 = inset_axes(axs[0,0], width="15%", height="15%", loc='lower left',bbox_transform=axs[0,0].transAxes,bbox_to_anchor=(0.1,0.25,1,1))
vax1 = inset_axes(axs[1,0], width="15%", height="15%", loc='lower left',bbox_transform=axs[1,0].transAxes,bbox_to_anchor=(0.1,0.25,1,1))

for vax in [vax0,vax1]:
    vax.plot(t,v,color=tab20c(4))
    vax.set_ylabel('Current',fontsize=6)
    vax.set_ylim(-3,3)
    vax.set_xlabel('Time',fontsize=6)

# extra formatting of insets
for pvax in [pax0,pax1,vax0,vax1]:
    pvax.set_xticks([])
    pvax.set_yticks([0])
    pvax.tick_params(axis='both', which='major', labelsize=6)
    pvax.axhline(0,linestyle='dashed',color='k')

# "sample" at +5m
samp_loc1 = 2.5
xi1 = np.argmin(np.abs(x-samp_loc1))
samp_loc2 = 5
xi2 = np.argmin(np.abs(x-samp_loc2))

# "sample" every 1 second
samp_times = np.arange(1,t.max(),1)
ti_list = [np.argmin(np.abs(t-samp_time)) for samp_time in samp_times]

mark1 = 'x'
mark2 = 'o'
col1 = tab20c(0)
col2 = tab20c(8)

for ti in ti_list:
    axs[0,0].plot(tt[ti,xi1],xx[ti,xi1],marker=mark1,mec=col1,mfc='None')
    axs[0,1].plot(t[ti],k_1noise_v[ti,xi1],marker=mark1,mec=col1,mfc='None')
    axs[1,0].plot(tt[ti,xi2],xx[ti,xi2],marker=mark2,mec=col2,mfc='None')
    axs[0,1].plot(t[ti],k_2noise_v[ti,xi2],marker=mark2,mec=col2,mfc='None')

axs[0,1].set_xlabel('Time (s)')
axs[0,1].set_ylabel('eDNA conc')

plt.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.95,hspace=0.3,wspace=0.2)
plt.show(block=False)
plt.pause(0.1)
# outfn = home + 'plots/displacement_kernel/figure_{:03d}.png'.format(tt)
# plt.savefig(outfn)
# plt.close()
