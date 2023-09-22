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
ylim = [0,0.55]

#define spatial grid
x = np.linspace(-10,10,101)
y = np.linspace(-10,10,101)

#define some coefficients
D = 1
alpha = 1
sigma = np.sqrt(2*D/alpha)
b = np.sqrt(D/alpha)


# solve for steady state analytical displacement kernel
c0 = 1/b
# c0 = np.sqrt(alpha/D)
# c1 = np.sqrt(alpha/D)
c1 = 1/(np.pi*(sigma**2))
c2 = 2/(sigma**2)

k1D = 0.5*c0 * np.exp(-c0*np.abs(x))

# k2D = np.zeros((len(y),len(x)))
# for i in range(len(x)):
#     for j in range(len(y)):
#         k2D[j,i] = 0.5*c0 * np.exp(-c0*np.sqrt(x[i]**2+y[j]**2))
#         # arg = np.sqrt(c2*(x[i]**2+y[j]**2))
#         # k2D[j,i] = c1*special.kn(0,arg)

# now solve for time dependent spin up
# initialize time grid
t = np.linspace(0,5,1000)

#solve for dx, dt
dx = x[1]-x[0]
# dy = y[1]-y[0]
dt = t[1]-t[0]

# initialize k solution
k1Dt = np.zeros((len(t),len(x)))
k2Dt = np.zeros((len(t),len(y),len(x)))

production = 10*np.max(k1D)

#condition at t=0 (index 50 when x=np.linspace(-50,50,101))
production_1D = np.zeros(len(x))
production_1D[50] = production
k1Dt[0,:] += production_1D*dt

# production_2D = np.zeros((len(y),len(x)))
# production_2D[50,50] = production#*1.5*np.pi
# k2Dt[0,:] += production_2D*dt

for tt in range(1,len(t)):
    #1D
    dkdx_up_1D = (k1Dt[tt-1,0:-2]-k1Dt[tt-1,1:-1])/dx
    dkdx_dn_1D = (k1Dt[tt-1,1:-1]-k1Dt[tt-1,2:])/dx
    diffusion_1D = D*(dkdx_up_1D-dkdx_dn_1D)/dx
    decay_1D = -k1Dt[tt-1,:]*alpha
    k1Dt[tt,1:-1] = k1Dt[tt-1,1:-1] + dt*(diffusion_1D + decay_1D[1:-1] + production_1D[1:-1])
    
    # #2D
    # dkdx_up_2D = (k2Dt[tt-1,:,0:-2]-k2Dt[tt-1,:,1:-1])/dx
    # dkdx_dn_2D = (k2Dt[tt-1,:,1:-1]-k2Dt[tt-1,:,2:])/dx
    # dkdy_up_2D = (k2Dt[tt-1,0:-2,:]-k2Dt[tt-1,1:-1,:])/dy
    # dkdy_dn_2D = (k2Dt[tt-1,1:-1,:]-k2Dt[tt-1,2:,:])/dy
    # diffusion_x_2D = 0.25*D*(dkdx_up_2D-dkdx_dn_2D)/dx
    # diffusion_y_2D = 0.25*D*(dkdy_up_2D-dkdy_dn_2D)/dy
    # decay_2D = -k2Dt[tt-1,:,:]*alpha*0
    # k2Dt[tt,1:-1,1:-1] = k2Dt[tt-1,1:-1,1:-1] + dt*(diffusion_x_2D[1:-1,:] + diffusion_y_2D[:,1:-1] + decay_2D[1:-1,1:-1] + production_2D[1:-1,1:-1])


#how many to plot
ntimes = 48
cmap = plt.get_cmap(cmo.cm.dense)
step = int(len(t)/ntimes)

# for time in range(0,ntimes):
for time in [ntimes-1]:
    fig,axs = plt.subplots(figsize=(10,8),nrows=2,ncols=2)
    axs[0,0].plot(x,k1D,color=tab20c(1),linestyle='solid')
    # axs[0,0].plot(x,k2D[50,:],color=tab20c(0),linestyle='dashed')
    axs[0,0].set_xlabel('x')
    axs[0,0].set_ylabel('k')
    axs[0,0].set_title('Analytical 1D displacement kernel w/ D=1, v=0')
    # axs[0,0].text(-8,0.8,'1D = solid\n2D = dashed'.format(tt),ha='left',va='center')
    axs[0,0].text(-8,0.8,'1D = solid'.format(tt),ha='left',va='center')
    # ylim = axs[0,0].get_ylim()


    # xx,yy = np.meshgrid(x,y)
    # p=axs[0,1].pcolormesh(xx,yy,k2D,cmap=cmo.cm.matter,vmin=vmin,vmax=vmax)
    # c=axs[0,1].contour(xx,yy,k2D,levels=np.arange(0,8,0.1),colors='k',linewidths=0.2)
    # axs[0,1].axhline(y[50],color=tab20c(0),linestyle='dashed')
    # axs[0,1].set_xlabel('x')
    # axs[0,1].set_ylabel('y')
    # axs[0,1].set_title('Analytical 2D displacement kernel w/ D=1, (u,v)=0')
    # cbaxes = inset_axes(axs[0,1], width="4%", height="60%", loc='lower right',bbox_transform=axs[0,1].transAxes,bbox_to_anchor=(0.1,-0.18,1,1))
    # cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
    # cb.set_label('k')
    # cb.add_lines(c)

    #time varying plots

    tt = int(time*step)
    # ct = int(time*256/ntimes)
    
    axs[1,0].plot(x,k1Dt[tt,:],color=tab20c(1),linestyle='solid')
    # axs[1,0].plot(x,k2Dt[tt,50,:],color=tab20c(0),linestyle='dashed')
    axs[1,0].set_xlabel('x')
    axs[1,0].set_ylabel('k')
    axs[1,0].set_title('Numercial 1D displacement kernel w/ D=1, v=0')
    axs[1,0].text(-8,0.4,'1D = solid\n'+r'time step = {:}'.format(tt),ha='left',va='center')
    # axs[1,0].set_ylim(ylim)
    
    # axs[1,1].pcolormesh(xx,yy,k2Dt[tt,:,:],cmap=cmo.cm.matter,vmin=vmin,vmax=vmax)
    # axs[1,1].contour(xx,yy,k2Dt[tt,:,:],levels=np.arange(0,8,0.1),colors='k',linewidths=0.2)
    # axs[1,1].axhline(y[50],color=tab20c(0),linestyle='dashed')
    # axs[1,1].set_xlabel('x')
    # axs[1,1].set_ylabel('y')
    # axs[1,1].set_title('Numerical 2D displacement kernel w/ D=1, (u,v)=0')

    for ax in axs[:,1]:
        ax.set_aspect(1)
    for ax in axs[:,0]:
        ax.set_ylim(ylim)
        
    # axs[0,2].plot(x,k2D[])

    plt.subplots_adjust(left=0.07,right=0.9,bottom=0.1,top=0.95,hspace=0.3,wspace=0.2)
    plt.show(block=False)
    plt.pause(0.1)
    # outfn = home + 'plots/displacement_kernel/figure_{:03d}.png'.format(tt)
    # plt.savefig(outfn)
    # plt.close()
