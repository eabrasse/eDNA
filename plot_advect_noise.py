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
k_constant = np.zeros((len(t),len(x)))
k_constant_v = np.zeros((len(t),len(x)))
k_pulse = np.zeros((len(t),len(x)))
k_pulse_v = np.zeros((len(t),len(x)))
k_noise = np.zeros((len(t),len(x)))
k_noise_v = np.zeros((len(t),len(x)))

#production scale
prod = 10

#condition at t=0 (index 50 when x=np.linspace(-50,50,101))
production_const = np.zeros((len(t),len(x)))
production_const[:,50] = prod

production_pulse = np.zeros((len(t),len(x)))
production_pulse[0,50] = prod*(t[-1]-t[0])

production_noise = np.zeros((len(t),len(x)))
# production_noise[:,50] = np.random.normal(prod,0.5*prod,len(t))
production_noise[:,50] = prod+prod*(np.cos(t*7*np.pi/(t[-1]-t[0]))+np.sin(t*2*np.pi/(t[-1]-t[0])))
production_noise[production_noise<0]=0
# norm = np.sum(production_noise)
# production_noise = production_noise

k_constant[0,:] += production_const[0,:]*dt
k_constant_v[0,:] += production_const[0,:]*dt

k_pulse[0,:] += production_pulse[0,:]
k_pulse_v[0,:] += production_pulse[0,:]

k_noise[0,:] += production_noise[0,:]*dt
k_noise_v[0,:] += production_noise[0,:]*dt

for tt in range(1,len(t)):
    #constant shedding, no current
    dkdx_up_const = (k_constant[tt-1,0:-2]-k_constant[tt-1,1:-1])/dx
    dkdx_dn_const = (k_constant[tt-1,1:-1]-k_constant[tt-1,2:])/dx
    diffusion_const = D*(dkdx_up_const-dkdx_dn_const)/dx
    decay_const = -k_constant[tt-1,:]*alpha
    k_constant[tt,1:-1] = k_constant[tt-1,1:-1] + dt*(diffusion_const + decay_const[1:-1] + production_const[tt,1:-1])
    
    #constant shedding, current
    dkdx_up_const_v = (k_constant_v[tt-1,0:-2]-k_constant_v[tt-1,1:-1])/dx
    dkdx_dn_const_v = (k_constant_v[tt-1,1:-1]-k_constant_v[tt-1,2:])/dx
    if v[tt-1]>=0:
        adv_const_v = v[tt-1]*dkdx_up_const_v
    elif v[tt-1]<0:
        adv_const_v = v[tt-1]*dkdx_dn_const_v
    diffusion_const_v = D*(dkdx_up_const_v-dkdx_dn_const_v)/dx
    decay_const_v = -k_constant_v[tt-1,:]*alpha
    k_constant_v[tt,1:-1] = k_constant_v[tt-1,1:-1] + dt*(adv_const_v + diffusion_const_v + decay_const_v[1:-1] + production_const[tt,1:-1])
    
    # pulse, no current
    dkdx_up_pulse = (k_pulse[tt-1,0:-2]-k_pulse[tt-1,1:-1])/dx
    dkdx_dn_pulse = (k_pulse[tt-1,1:-1]-k_pulse[tt-1,2:])/dx
    diffusion_pulse = D*(dkdx_up_pulse-dkdx_dn_pulse)/dx
    decay_pulse = -k_pulse[tt-1,:]*alpha
    k_pulse[tt,1:-1] = k_pulse[tt-1,1:-1] + dt*(diffusion_pulse + decay_pulse[1:-1]+ production_pulse[tt,1:-1])
    
    # pulse, current
    dkdx_up_pulse_v = (k_pulse_v[tt-1,0:-2]-k_pulse_v[tt-1,1:-1])/dx
    dkdx_dn_pulse_v = (k_pulse_v[tt-1,1:-1]-k_pulse_v[tt-1,2:])/dx
    if v[tt-1]>=0:
        adv_pulse_v = v[tt-1]*dkdx_up_pulse_v
    elif v[tt-1]<0:
        adv_pulse_v = v[tt-1]*dkdx_dn_pulse_v
    diffusion_pulse_v = D*(dkdx_up_pulse_v-dkdx_dn_pulse_v)/dx
    decay_pulse_v = -k_pulse_v[tt-1,:]*alpha
    k_pulse_v[tt,1:-1] = k_pulse_v[tt-1,1:-1] + dt*(adv_pulse_v + diffusion_pulse_v + decay_pulse_v[1:-1]+ production_pulse[tt,1:-1])
    
    # noise, no current
    dkdx_up_noise = (k_noise[tt-1,0:-2]-k_noise[tt-1,1:-1])/dx
    dkdx_dn_noise = (k_noise[tt-1,1:-1]-k_noise[tt-1,2:])/dx
    diffusion_noise = D*(dkdx_up_noise-dkdx_dn_noise)/dx
    decay_noise = -k_noise[tt-1,:]*alpha
    k_noise[tt,1:-1] = k_noise[tt-1,1:-1] + dt*(diffusion_noise + decay_noise[1:-1]+ production_noise[tt,1:-1])
    
    # pulse, current
    dkdx_up_noise_v = (k_noise_v[tt-1,0:-2]-k_noise_v[tt-1,1:-1])/dx
    dkdx_dn_noise_v = (k_noise_v[tt-1,1:-1]-k_noise_v[tt-1,2:])/dx
    if v[tt-1]>=0:
        adv_noise_v = v[tt-1]*dkdx_up_noise_v
    elif v[tt-1]<0:
        adv_noise_v = v[tt-1]*dkdx_dn_noise_v
    diffusion_noise_v = D*(dkdx_up_noise_v-dkdx_dn_noise_v)/dx
    decay_noise_v = -k_noise_v[tt-1,:]*alpha
    k_noise_v[tt,1:-1] = k_noise_v[tt-1,1:-1] + dt*(adv_noise_v + diffusion_noise_v + decay_noise_v[1:-1]+ production_noise[tt,1:-1])
    

#plotting params
cmap = plt.get_cmap(cmo.cm.dense)
xx,tt = np.meshgrid(x,t)

#generate figure
fig,axs = plt.subplots(figsize=(11,8),nrows=2,ncols=3)

for ax,k in [[axs[0,0],k_constant],[axs[1,0],k_constant_v],[axs[0,1],k_pulse],[axs[1,1],k_pulse_v],[axs[0,2],k_noise],[axs[1,2],k_noise_v]]:
    p=ax.pcolormesh(tt,xx,k,cmap=cmo.cm.matter,vmin=vmin,vmax=vmax)
    ax.axhline(0,linestyle='dashed',color='k')
    ax.set_aspect(0.25)
for ax in axs[1,:]:
    ax.set_xlabel('Time (s)')
for ax in axs[:,0]:
    ax.set_ylabel('Distance (m)')

axs[0,0].text(0.1,0.8,'Constant eDNA shedding\nNo current',ha='left',va='center',transform=axs[0,0].transAxes)
axs[1,0].text(0.1,0.8,'Constant eDNA shedding\nCurrent',ha='left',va='center',transform=axs[1,0].transAxes)
axs[0,1].text(0.1,0.8,'eDNA pulse\nNo current',ha='left',va='center',transform=axs[0,1].transAxes)
axs[1,1].text(0.1,0.8,'eDNA pulse\nCurrent',ha='left',va='center',transform=axs[1,1].transAxes)
axs[0,2].text(0.1,0.8,'eDNA noise\nNo current',ha='left',va='center',transform=axs[0,2].transAxes)
axs[1,2].text(0.1,0.8,'eDNA noise\nCurrent',ha='left',va='center',transform=axs[1,2].transAxes)

cbaxes = inset_axes(axs[0,-1], width="4%", height="60%", loc='center right',bbox_transform=axs[0,-1].transAxes,bbox_to_anchor=(0.1,0.,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')

pax0 = inset_axes(axs[0,0], width="15%", height="15%", loc='upper right',bbox_transform=axs[0,0].transAxes,bbox_to_anchor=(0.,-0.05,1,1))
pax1 = inset_axes(axs[1,0], width="15%", height="15%", loc='upper right',bbox_transform=axs[1,0].transAxes,bbox_to_anchor=(0.,-0.05,1,1))
pax2 = inset_axes(axs[0,1], width="15%", height="15%", loc='upper right',bbox_transform=axs[0,1].transAxes,bbox_to_anchor=(0.,-0.05,1,1))
pax3 = inset_axes(axs[1,1], width="15%", height="15%", loc='upper right',bbox_transform=axs[1,1].transAxes,bbox_to_anchor=(0.,-0.05,1,1))
pax4 = inset_axes(axs[0,2], width="15%", height="15%", loc='upper right',bbox_transform=axs[0,2].transAxes,bbox_to_anchor=(0.,-0.05,1,1))
pax5 = inset_axes(axs[1,2], width="15%", height="15%", loc='upper right',bbox_transform=axs[1,2].transAxes,bbox_to_anchor=(0.,-0.05,1,1))

for pax in [pax0,pax1]:
    pax.plot(t,production_const[:,50],color=tab20c(0))
    pax.set_ylim(-5,20)
for pax in [pax2,pax3]:
    pax.plot(t,production_pulse[:,50],color=tab20c(0))
    pax.set_ylim([-15,60])
for pax in [pax4,pax5]:
    pax.plot(t,production_noise[:,50],color=tab20c(0))
    pax.set_ylim([-10,30])
for pax in [pax0,pax1,pax2,pax3,pax4,pax5]:
    pax.set_ylabel('eDNA production',fontsize=6)

vax0 = inset_axes(axs[0,0], width="15%", height="15%", loc='upper right',bbox_transform=axs[0,0].transAxes,bbox_to_anchor=(0.,-0.25,1,1))
vax1 = inset_axes(axs[1,0], width="15%", height="15%", loc='upper right',bbox_transform=axs[1,0].transAxes,bbox_to_anchor=(0.,-0.25,1,1))
vax2 = inset_axes(axs[0,1], width="15%", height="15%", loc='upper right',bbox_transform=axs[0,1].transAxes,bbox_to_anchor=(0.,-0.25,1,1))
vax3 = inset_axes(axs[1,1], width="15%", height="15%", loc='upper right',bbox_transform=axs[1,1].transAxes,bbox_to_anchor=(0.,-0.25,1,1))
vax4 = inset_axes(axs[0,2], width="15%", height="15%", loc='upper right',bbox_transform=axs[0,2].transAxes,bbox_to_anchor=(0.,-0.25,1,1))
vax5 = inset_axes(axs[1,2], width="15%", height="15%", loc='upper right',bbox_transform=axs[1,2].transAxes,bbox_to_anchor=(0.,-0.25,1,1))

for vax in [vax0,vax2,vax4]:
    vax.plot(t,v0,color=tab20c(4))
for vax in [vax1,vax3,vax5]:
    vax.plot(t,v,color=tab20c(4))
for vax in [vax0,vax1,vax2,vax3,vax4,vax5]:
    vax.set_ylabel('Current',fontsize=6)
    vax.set_ylim(-3,3)
    vax.set_xlabel('Time',fontsize=6)

for pvax in [pax0,pax1,pax2,pax3,pax4,pax5,vax0,vax1,vax2,vax3,vax4,vax5]:
    pvax.set_xticks([])
    pvax.set_yticks([0])
    pvax.tick_params(axis='both', which='major', labelsize=6)
    pvax.axhline(0,linestyle='dashed',color='k')
    

plt.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.95,hspace=0.3,wspace=0.2)
plt.show(block=False)
plt.pause(0.1)
# outfn = home + 'plots/displacement_kernel/figure_{:03d}.png'.format(tt)
# plt.savefig(outfn)
# plt.close()
