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
vmax = 0.8
ylim = [0,0.55]

#define spatial grid
x = np.linspace(-10,10,101)

#define some coefficients
D = 1
alpha = 1
sigma = np.sqrt(2*D/alpha)
b = np.sqrt(D/alpha)

# now solve for time dependent spin up
# initialize time grid
t = np.linspace(0,5,1000)

#solve for dx, dt
dx = x[1]-x[0]
dt = t[1]-t[0]

# initialize k solution
k_constant = np.zeros((len(t),len(x)))
k_pulse = np.zeros((len(t),len(x)))

production = 10

#condition at t=0 (index 50 when x=np.linspace(-50,50,101))
production_vec = np.zeros(len(x))
production_vec[50] = production
k_constant[0,:] += production_vec*dt
k_pulse[0,:] += production_vec*(t[-1]-t[0])

for tt in range(1,len(t)):
    #1D
    dkdx_up_const = (k_constant[tt-1,0:-2]-k_constant[tt-1,1:-1])/dx
    dkdx_dn_const = (k_constant[tt-1,1:-1]-k_constant[tt-1,2:])/dx
    diffusion_const = D*(dkdx_up_const-dkdx_dn_const)/dx
    decay_const = -k_constant[tt-1,:]*alpha
    k_constant[tt,1:-1] = k_constant[tt-1,1:-1] + dt*(diffusion_const + decay_const[1:-1] + production_vec[1:-1])
    
    dkdx_up_pulse = (k_pulse[tt-1,0:-2]-k_pulse[tt-1,1:-1])/dx
    dkdx_dn_pulse = (k_pulse[tt-1,1:-1]-k_pulse[tt-1,2:])/dx
    diffusion_pulse = D*(dkdx_up_pulse-dkdx_dn_pulse)/dx
    decay_pulse = -k_pulse[tt-1,:]*alpha
    k_pulse[tt,1:-1] = k_pulse[tt-1,1:-1] + dt*(diffusion_pulse + decay_pulse[1:-1])
    

#plotting params
cmap = plt.get_cmap(cmo.cm.dense)
xx,tt = np.meshgrid(x,t)

#generate figure
fig,axs = plt.subplots(figsize=(10,5),nrows=1,ncols=2)
p=axs[0].pcolormesh(tt,xx,k_constant,cmap=cmo.cm.matter,vmin=vmin,vmax=vmax)
axs[0].set_xlabel('Time (s)')
axs[0].set_ylabel('Distance (m)')
axs[0].axhline(0,linestyle='dashed',color='k')
axs[0].text(0.1,0.8,'Constant eDNA shedding',ha='left',va='center',transform=axs[0].transAxes)

prod_const = production*dt*np.ones(t.shape)
p_axes0 = inset_axes(axs[0], width="15%", height="15%", loc='upper right',bbox_transform=axs[0].transAxes,bbox_to_anchor=(0.,-0.05,1,1))
p_axes0.plot(t,prod_const,color=tab20c(0))
p_axes0.set_xticks([])
p_axes0.set_yticks([0])
p_axes0.tick_params(axis='both', which='major', labelsize=6)
p_axes0.set_xlabel('Time',fontsize=6)
p_axes0.set_ylabel('eDNA production',fontsize=6)
p_axes0.set_ylim(-0.05,0.2)
p_axes0.axhline(0,linestyle='dashed',color='k')

axs[1].pcolormesh(tt,xx,k_pulse,cmap=cmo.cm.matter,vmin=vmin,vmax=vmax)
axs[1].set_xlabel('Time (s)')
# axs[1].set_ylabel('Distance')
axs[1].axhline(0,linestyle='dashed',color='k')
axs[1].text(0.1,0.8,'eDNA pulse',ha='left',va='center',transform=axs[1].transAxes)

prod_pulse = np.zeros(t.shape)
prod_pulse[0] = production*(t[-1]-t[0])
p_axes1 = inset_axes(axs[1], width="15%", height="15%", loc='upper right',bbox_transform=axs[1].transAxes,bbox_to_anchor=(0.,-0.05,1,1))
p_axes1.plot(t,prod_pulse,color=tab20c(0))
p_axes1.set_xticks([])
p_axes1.set_yticks([0])
p_axes1.tick_params(axis='both', which='major', labelsize=6)
p_axes1.set_xlabel('Time',fontsize=6)
p_axes1.set_ylabel('eDNA production',fontsize=6)
p_axes1.set_ylim([-15,60])
p_axes1.axhline(0,linestyle='dashed',color='k')

cbaxes = inset_axes(axs[1], width="4%", height="60%", loc='center right',bbox_transform=axs[1].transAxes,bbox_to_anchor=(0.1,0.,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')

for ax in axs:
    ax.set_aspect(0.25)

plt.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.95,hspace=0.3,wspace=0.2)
plt.show(block=False)
plt.pause(0.1)
# outfn = home + 'plots/displacement_kernel/figure_{:03d}.png'.format(tt)
# plt.savefig(outfn)
# plt.close()
