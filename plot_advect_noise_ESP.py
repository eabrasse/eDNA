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
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import pytz
import pickle
import calendar
from matplotlib.gridspec import GridSpec

def toTimestamp(d):
  return calendar.timegm(d.timetuple())

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'
tab20c = plt.get_cmap('tab20c',20)
tz = pytz.utc

#load in tidal velocity
vel_fn = home+'data/HC_dolph_tidal_current_ATTideGauge_9445133.p'
D = pickle.load(open(vel_fn,'rb'))
dt_list = D['dt_list'][:]
dt_list = [tz.localize(dt) for dt in dt_list]
dt_list2 = np.array([toTimestamp(dt) for dt in dt_list])
u0 = D['u'][:]
v0 = D['v'][:]
# v = np.sqrt(u0**2+v0**2)
w0 = u0 + 1j*v0
w = w0-w0.mean()
cov = np.cov(np.imag(w),np.real(w))
term1=cov[0,0]+cov[1,1]
term2=np.sqrt((cov[0,0]-cov[1,1])**2 + 4*cov[1,0]**2)
major=np.sqrt(.5*(term1+term2))
minor=np.sqrt(.5*(term1-term2))
theta = 0.5*np.arctan2(2*cov[1,0],(cov[0,0]-cov[1,1]))
w_pa = w0*np.exp(1j*theta)
v=np.imag(w_pa)

# load in samples data
ESP_fn = home+'data/03_ESP1_Feb2023_hourly.csv'
df = pd.read_csv(ESP_fn,sep=',',engine='python')
df['datetime']=pd.to_datetime(df.date)
# dt_DNA = [td.seconds for td in np.diff(df.dropna().datetime)]
# dDNA = np.diff(df.dropna().PB_quantity_mean)
# LHS = dDNA/dt_DNA

times = pd.to_datetime(df.dropna().datetime.values)
times = [tz.localize(time) for time in times]
times = [toTimestamp(time) for time in times]
tmin = min(times)
tmax = max(times)
# times = [times[i] + 0.5*dt_DNA[i] for i in range(len(dt_DNA))]
velmask = np.array([(dt>tmin)&(dt<tmax) for dt in dt_list2])
dt_list = [dt_list[i] for i in range(len(dt_list)) if velmask[i]]
dt_list2 = dt_list2[velmask]
v = v[velmask]
production = np.interp(dt_list2,np.array(times),df.dropna().PB_quantity_mean)


ylim = [0,0.55]

# #define spatial grid
x = np.linspace(-100000,100000,51)

#define some coefficients
D = 0
alpha = 5.5e-6
sigma = np.sqrt(2*D/alpha)
b = np.sqrt(D/alpha)

# now solve for time dependent spin up

#solve for dx, dt
dx = x[1]-x[0]
dt = dt_list2[1]-dt_list2[0]

# initialize k solution
k = np.zeros((len(dt_list),len(x)))
adv = np.zeros((len(dt_list),len(x)))

k[0,25] = production[0]
for tt in range(1,len(dt_list)):
    
    # noise, current
    dkdx_up = (k[tt-1,0:-2]-k[tt-1,1:-1])/dx
    dkdx_dn = (k[tt-1,1:-1]-k[tt-1,2:])/dx
    if v[tt-1]>=0:
        adv[tt,1:-1] = v[tt-1]*dkdx_up
    elif v[tt-1]<0:
        adv[tt,1:-1] = v[tt-1]*dkdx_dn
    diffusion = D*(dkdx_up-dkdx_dn)/dx
    decay= -k[tt-1,:]*alpha
    k[tt,1:-1] = k[tt-1,1:-1] + dt*(adv[tt,1:-1] + diffusion + decay[1:-1])
    k[tt,25] += dt*production[tt]
    

#plotting params
cmap = plt.get_cmap(cmo.cm.dense)
xx,tt = np.meshgrid(x,dt_list2)

#generate figure
fig = plt.figure(figsize=(12,8))
gs = GridSpec(3,3)
ax = fig.add_subplot(gs[:,:2])
p=ax.pcolormesh(dt_list,0.001*x,k.T,cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e8),shading='nearest')
ax.axhline(0,linestyle='dashed',color='k')
# ax.set_aspect(0.25)
ax.set_xlabel('Date')
ax.set_ylabel('Distance (km)')
cbaxes = inset_axes(ax, width="40%", height="6%", loc='lower left',bbox_transform=ax.transAxes,bbox_to_anchor=(0.1,0.1,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='horizontal')

pax = inset_axes(ax, width="15%", height="15%", loc='upper right',bbox_transform=ax.transAxes,bbox_to_anchor=(0.,-0.05,1,1))
pax.plot(dt_list,production,color=tab20c(0))
pax.set_yscale('log')
pax.set_ylabel('eDNA production',fontsize=6)

vax = inset_axes(ax, width="15%", height="15%", loc='upper right',bbox_transform=ax.transAxes,bbox_to_anchor=(0.,-0.25,1,1))
vax.plot(dt_list,v,color=tab20c(4))
vax.set_ylabel('Current',fontsize=6)
vax.set_ylim([-0.5,0.5])
vax.set_xlabel('Time',fontsize=6)
vax.set_yticks([0])
vax.axhline(0,linestyle='dashed',color='k')

for pvax in [pax,vax]:
    pvax.set_xticks([])
    pvax.tick_params(axis='both', which='major', labelsize=6)
ax.text(0.1,0.9,'a) Modeled DNA distribution',transform=ax.transAxes)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d %H:%m"))
plt.setp( ax.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')

t0=5
ax_ts = fig.add_subplot(gs[0,2])
ax_ts.scatter(dt_list[t0:],adv[t0:,25],c=times,cmap='rainbow')
ax_ts.set_xlabel('Date')
ax_ts.set_ylabel('copies/s')
ax_ts.text(0.1,0.9,'b) Advection',transform=ax_ts.transAxes)
ax_ts.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d %H:%m"))
plt.setp( ax_ts.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')

advlim = -5000
ax_av = fig.add_subplot(gs[1,2])
ax_av.scatter(v[t0:],adv[t0:,25],c=times,cmap='rainbow')
ax_av.set_ylabel('Advection (copies/s)')
ax_av.set_xlabel('Current (m/s)')
ax_av.text(0.1,0.9,'c) Advection vs Velocity',transform=ax_av.transAxes)
p = np.polyfit(v[t0:][adv[t0:,25]>advlim],adv[t0:,25][adv[t0:,25]>advlim],1)
px = v[t0:][adv[t0:,25]>advlim]
px.sort()
py = p[0]*px + p[1]
ax_av.plot(px,py,linestyle='dashed',color='k')
r = np.corrcoef(adv[t0:,25][adv[t0:,25]>advlim],v[t0:][adv[t0:,25]>advlim])[0,1]
Rsquared = r**2
ax_av.text(0.9,0.1,r'$\mathrm{R}^{2}$='+f'{Rsquared:.2}',transform=ax_av.transAxes,ha='right')

dDNA = k[1:,25]-k[:-1,25]
dr = v[:-1]*dt
dDNAdr = dDNA/dr
vdDNAdr = dDNAdr*v[1:]
ax_ak = fig.add_subplot(gs[2,2])
ax_ak.scatter(vdDNAdr[(t0-1):],adv[t0:,25],c=times,cmap='rainbow')
ax_ak.set_ylabel('Advection (copies/s)')
ax_ak.set_xlabel(r'$v\frac{\partial DNA}{\partial r}$')
ax_ak.text(0.1,0.9,'d) Advection vs estimated gradient times current',transform=ax_ak.transAxes)
p = np.polyfit(vdDNAdr[(t0-1):][adv[t0:,25]>advlim],adv[t0:,25][adv[t0:,25]>advlim],1)
px = vdDNAdr[(t0-1):][adv[t0:,25]>advlim]
px.sort()
py = p[0]*px + p[1]
ax_ak.plot(px,py,linestyle='dashed',color='k')
r = np.corrcoef(adv[t0:,25][adv[t0:,25]>advlim],vdDNAdr[(t0-1):][adv[t0:,25]>advlim])[0,1]
Rsquared = r**2
ax_ak.text(0.9,0.1,r'$\mathrm{R}^{2}$='+f'{Rsquared:.2}',transform=ax_ak.transAxes,ha='right')

# for ax in ax_av,ax_ak:
#     ax.set_xlim([advlim,1000])

plt.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.95,hspace=0.4,wspace=0.4)
plt.show(block=False)
plt.pause(0.1)
# outfn = home + 'plots/displacement_kernel/figure_{:03d}.png'.format(tt)
# plt.savefig(outfn)
# plt.close()
