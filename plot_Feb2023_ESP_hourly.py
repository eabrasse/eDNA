#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Just testing to see if the hourly compression worked
"""

# importing modules, the beginning of all python code
import os
import sys
import matplotlib
matplotlib.use('macosx')
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

def toTimestamp(d):
  return calendar.timegm(d.timetuple())

# some helpful plotting commands
plt.close('all')
tab10 = plt.get_cmap('tab10',10)
rainbow = plt.get_cmap('rainbow')
props = dict(boxstyle='round', fc='white',ec='None',alpha=.5)

# load in data
home = '/Users/elizabethbrasseale/Projects/eDNA/data/'
data_fn = home+'03_ESP1_Feb2023_hourly.csv'

# generate figure and axis handles
# fig,axs = plt.subplots(nrows=3, ncols=1,figsize=(8,6))
fig = plt.figure(figsize=(7,10))
gs = GridSpec(6,1)
ax1 = fig.add_subplot(gs[2:-2,0])

df = pd.read_csv(data_fn,sep=',',engine='python')
df['datetime']=pd.to_datetime(df.date)
# df.datetime = df.datetime.dt.tz_localize(tz='UTC')


ax1.scatter(df.datetime,df.PB_quantity_mean,c=tab10(0))#edgecolors=tab10(0),c='None')
ax1.set_yscale('log')
ax1.grid()
# ax1.set_xlabel('Time (UTC)')
ax1.set_ylabel(r'Conc (copies/$\mathrm{\mu L}$)')
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d %H:%m"))
# plt.setp( ax1.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax1.text(0.1,0.9,'b) ESP samples near net pen',transform=ax1.transAxes,bbox=props)
# Tues2pm = datetime(2023,1,31,14,0,0)
# pdt = pytz.timezone('America/Vancouver')
# Tues2pm = pdt.localize(Tues2pm)
# Tues2pm = Tues2pm.astimezone(pytz.utc)
# ax.axvline(x=Tues2pm,color=tab10(1))
# ax.text(Tues2pm,1e4,'Tues 2pm',color=tab10(1))
tz = pytz.utc
xlim = [tz.localize(datetime(2023,1,30,23)),tz.localize(datetime(2023,2,2,2))]


# add dolphin count
dolphcount_fn = home+'dolphin_presence.csv'
dcf = pd.read_csv(dolphcount_fn,sep=',',engine='python')
dcf['datetime']=pd.to_datetime(dcf.date+' '+dcf.time)
dcf= dcf[dcf.ndolphin_net.notnull()]
dcf.datetime=dcf.datetime.dt.tz_localize(tz='America/Vancouver')
dcf.datetime=dcf.datetime.dt.tz_convert('UTC')


ax2 = fig.add_subplot(gs[-1,0])
ax2.bar(dcf.datetime+timedelta(minutes=30),dcf.ndolphin_net,color=tab10(1),width=0.03)
ax2.set_ylabel('# dolphins')
ax2.set_yticks([0,1,2,3])
ax2.text(0.1,0.9,'d) Number of dolphins in net pen',transform=ax2.transAxes,bbox=props)
# ax2.set_xlim(xlim)
ax2.grid()

# add tides
mask_fn = home+'HC_dolph_tidal_current_ATTideGauge_ks0201.p'
D = pickle.load(open(mask_fn,'rb'))

dt_list = D['dt_list'][:]
dt_list = [tz.localize(dt) for dt in dt_list]

u = D['u'][:]
v = D['v'][:]

tvel = np.sqrt(u**2+v**2)
# tvel = np.sqrt(u**2)
ax3 = fig.add_subplot(gs[-2,0])
ax3.plot(dt_list,tvel,c=tab10(2))
ax3.set_ylabel('velocity (m/s)')
ax3.text(0.1,0.9,'c) Tidal current speed',transform=ax3.transAxes,bbox=props)
# ax3.set_xlim(xlim)
ax3.grid()


# # add tide + dolphin number
# tmin = max([min(wd_dt_list),dcf.datetime.min()])
# tmax = min([max(wd_dt_list),dcf.datetime.max()])
# wd_dt_list_sub = [wd_dt for wd_dt in wd_dt_list if wd_dt>tmin and wd_dt<tmax]
# zeta_sub = [zeta[i] for i in range(len(wd_dt_list)) if wd_dt_list[i]>tmin and wd_dt_list[i]<tmax]
# dcf_sub = dcf[(dcf.datetime>tmin)&(dcf.datetime<tmax)]
# tidedolph = dcf_sub.ndolphin_net.values+zeta_sub

# ax4.plot(wd_dt_list_sub,tidedolph,c=tab10(3))
# ax4.set_ylabel('dolphin + meters (?)')
# ax4.grid()
# ax4.text(0.1,0.9,r'b) Tides + Dolphins',transform=ax4.transAxes,bbox=props)

ax4 = fig.add_subplot(gs[:2,0])

#calculate terms
dt_DNA = [td.seconds for td in np.diff(df.datetime)]
dDNA = np.diff(df.PB_quantity_mean)
LHS = dDNA/dt_DNA
# LHS = 0.5*[LHS[:-1]+LHS[1:]]
LHS_dt = df.datetime[:-1]+timedelta(minutes=30)

shedding_rate = 1e0 #guess
shedding0 = dcf.ndolphin_net*shedding_rate

# advection includes DNA & velocity which are currently on different time schedules
times = pd.to_datetime(df.dropna().datetime.values)
times = [tz.localize(time) for time in times]
times = [toTimestamp(time) for time in times]
dt_list2 = np.array([toTimestamp(dt) for dt in dt_list])
tmin = min(times)
tmax = max(times)
velmask = np.array([(dt>tmin)&(dt<tmax) for dt in dt_list2])
vel = np.interp(np.array(times),dt_list2[velmask],tvel[velmask])

dcf_times = pd.to_datetime(dcf.datetime.values)
dcf_times = [tz.localize(time) for time in dcf_times]
dcf_times = [toTimestamp(time) for time in dcf_times]
shedmask = np.array([(dt>tmin)&(dt<tmax) for dt in dcf_times])
shedding = np.interp(np.array(times),np.array(dcf_times)[shedmask],shedding0.values[shedmask])

# vel = []
# shedding = []
# # dt_list_sub = []
# for time in df.datetime:
#     td_list = [np.abs(time-dt) for dt in dt_list]
#     tind = np.argmin(td_list)
#     vel.append(tvel[tind])
#
#     td_list = [np.abs(time-dt) for dt in dcf.datetime]
#     tind = np.argmin(td_list)
#     shedding.append(shedding0.values[tind])
#     # dt_list_sub.append(dt_list[tind])
advection = -df.dropna().PB_quantity_mean.values*vel/1000

sinking_rate = 0 #guess
sinking = -df.dropna().PB_quantity_mean.values*sinking_rate

# decay_rate = np.array([0.02*dt_DNA0/3600 for dt_DNA0 in dt_DNA])
# decay_rate = 0.5*(decay_rate[1:]+decay_rate[:-1])
decay_rate = 0.02/3600
decay = -df.dropna().PB_quantity_mean.values[1:-1]*decay_rate

RHS = shedding[1:-1] + advection[1:-1] + decay

ax4.plot(LHS_dt,LHS,c=tab10(0),label=r'LHS ($\frac{\Delta\mathrm{DNA}}{\Delta t}$)')
ax4.plot(df.dropna().datetime.values,shedding,c=tab10(1),ls='dashed',label=r'shedding ($N_{dolph}R_{shed}$)')
ax4.plot(df.dropna().datetime.values,advection,c=tab10(2),ls='dashed',label=r'advection ($-V\frac{\mathrm{DNA}}{\Delta x}$)')
ax4.plot(df.dropna().datetime.values[1:-1],decay,c=tab10(3),ls='dashed',label=r'decay ($-k \mathrm{DNA}$)')
ax4.plot(df.dropna().datetime.values[1:-1],RHS,c=tab10(4),label=r'$\Sigma$ RHS')

ax4.set_ylabel(r'copies $\mathrm{second}^{-1}$')
ax4.grid()
ax4.text(0.1,0.9,r'a) Balance',transform=ax4.transAxes,bbox=props)
ax4.legend()
# ax4.set_yscale('log')

for ax in ax1,ax4,ax3:
    ax.set_xticklabels(['' for xt in ax.get_xticklabels()])
    ax.set_xlabel('')
ax2.set_xlabel('Time (UTC)')
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d %H:%m"))
plt.setp( ax2.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')

for ax in ax1,ax2,ax3,ax4:
    ax.set_xlim(xlim)

# show plot
# plt.tight_layout()
plt.subplots_adjust(top=0.98,hspace=0.1)
plt.show(block=False)
plt.pause(0.1)
