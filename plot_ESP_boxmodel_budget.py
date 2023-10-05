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
smdata_fn = home+'process_term_mean.csv'


# generate figure and axis handles
fig,axs = plt.subplots(nrows=3, ncols=1,figsize=(8,12))
# fig,axs = plt.subplots(3,1)
ax = axs[0]


df = pd.read_csv(data_fn,sep=',',engine='python')
df['datetime']=pd.to_datetime(df.date)
tz = pytz.utc
xlim = [tz.localize(datetime(2023,1,30,23)),tz.localize(datetime(2023,2,2,2))]
smdf = pd.read_csv(smdata_fn,sep=',',engine='python')
df['log_process_mean']=smdf['process_term'][:]
DNA_mean = np.mean(np.log(df.PB_quantity_mean))
df['process_mean'] = np.exp(df['log_process_mean']+DNA_mean)


# add dolphin count
dolphcount_fn = home+'dolphin_presence.csv'
dcf = pd.read_csv(dolphcount_fn,sep=',',engine='python')
dcf['datetime']=pd.to_datetime(dcf.date+' '+dcf.time)
dcf= dcf[dcf.ndolphin_net.notnull()]
dcf.datetime=dcf.datetime.dt.tz_localize(tz='America/Vancouver')
dcf.datetime=dcf.datetime.dt.tz_convert('UTC')


# add tides
mask_fn = home+'HC_dolph_tidal_current_ATTideGauge_ks0201.p'
D = pickle.load(open(mask_fn,'rb'))

dt_list = D['dt_list'][:]
dt_list = [tz.localize(dt) for dt in dt_list]

u = D['u'][:]
v = D['v'][:]

tvel = np.sqrt(u**2+v**2)


#calculate terms
dt_DNA = [td.seconds for td in np.diff(df.dropna().datetime)]
# dDNA = np.diff(df.dropna().PB_quantity_mean)
dDNA = np.diff(df.dropna().process_mean)
LHS = dDNA/dt_DNA
# LHS = 0.5*[LHS[:-1]+LHS[1:]]
LHS_dt = df.dropna().datetime[:-1]+timedelta(minutes=30)

# shedding_rate = 0.017*1e0 #guess
shedding_rate = 0.25
shedding = dcf.ndolphin_net*shedding_rate

# advection includes DNA & velocity which are currently on different time schedules
times = pd.to_datetime(df.dropna().datetime.values)
times = [tz.localize(time) for time in times]
times = [toTimestamp(time) for time in times]
dt_list2 = np.array([toTimestamp(dt) for dt in dt_list])
tmin = min(times)
tmax = max(times)
velmask = np.array([(dt>tmin)&(dt<tmax) for dt in dt_list2])
# vel = np.interp(np.array(times),dt_list2[velmask],tvel[velmask])

dcf_times = pd.to_datetime(dcf.datetime.values)
dcf_times = [tz.localize(time) for time in dcf_times]
dcf_times = [toTimestamp(time) for time in dcf_times]
shedmask = np.array([(dt>tmin)&(dt<tmax) for dt in dcf_times])
# shedding = np.interp(np.array(times),np.array(dcf_times)[shedmask],shedding0.values[shedmask])

V = 1
# PB_quantity_mean = V*np.interp(dt_list2[velmask],np.array(times),df.dropna().PB_quantity_mean.values)
process_term_mean = V*np.interp(dt_list2[velmask],np.array(times),df.dropna().process_mean.values)

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
# advection = -df.dropna().PB_quantity_mean.values*vel/1000
# advection = - PB_quantity_mean*tvel[velmask]/1000
C=1.3e-3
advection = - C * tvel[velmask] * process_term_mean

sinking_rate = 0 #guess
sinking = -process_term_mean*sinking_rate
# 
decay_rate = 1e-4
# decay_rate = 0.02/3600
# decay = -df.dropna().PB_quantity_mean.values[1:-1]*decay_rate
decay = -process_term_mean*decay_rate

RHS = shedding[shedmask] + advection + decay

ax.plot(LHS_dt,LHS,c=tab10(0),label=r'LHS ($\frac{\Delta\mathrm{DNA}}{\Delta t}$)',lw=2)
ax.plot(np.array(dt_list)[velmask],shedding[shedmask],c=tab10(1),ls='dashed',label=r'Shed ($N_{dolph}R_{shed}$)',alpha=0.5)
ax.plot(np.array(dt_list)[velmask],advection,c=tab10(2),ls='dashed',label=r'Adv ($-v_{r}\frac{\mathrm{DNA}}{\Delta r}$)',alpha=0.5)
# ax.plot(np.array(dt_list)[velmask],advection,c=tab10(2),ls='dashed',label=r'Adv ($-A v_{r}$)',alpha=0.5)
ax.plot(np.array(dt_list)[velmask],decay,c=tab10(3),ls='dashed',label=r'Decay ($-k \mathrm{DNA}$)',alpha=0.5)
ax.plot(np.array(dt_list)[velmask],RHS,c=tab10(4),label=r'RHS = Shed + Adv + Decay',lw=2)

ax.set_ylabel(r'copies $\mathrm{second}^{-1}$')
ax.grid()
ax.text(0.1,0.9,r'dDNA/dt Equation balance',transform=ax.transAxes,bbox=props)
ax.legend()
# ax.set_yscale('log')


ax.set_xlabel('Time (UTC)')
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d %-H:%m"))
plt.setp( ax.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax.set_xlim(xlim)


ax1 = axs[1]
DNA_0 = df.dropna().process_mean.values[0]
dt = np.array(dt_list)[velmask][1]-np.array(dt_list)[velmask][0]
RHS_cs = np.cumsum(np.insert(3600*RHS.values,0,DNA_0))
dt_list_cs = np.array(dt_list)[velmask]
dt_list_cs = np.insert(dt_list_cs,0,(dt_list_cs[0]-timedelta(hours=1)))
# print(DNA_0)
# RHS_cs+=df.process_mean[0]
# LHS_cs = df.process_mean
LHS_cs = np.cumsum(np.insert(dDNA,0,DNA_0))

ax1.plot(df.dropna().datetime,LHS_cs,c=tab10(0),label=r'LHS ($\frac{\Delta\mathrm{DNA}}{\Delta t}$)',lw=2)
ax1.plot(dt_list_cs,RHS_cs,c=tab10(4),label=r'RHS = Shed + Adv + Decay',lw=2)
ax1.set_ylabel(r'copies $\mathrm{second}^{-1}$')
ax1.grid()
ax1.text(0.1,0.9,r'DNA Equation balance',transform=ax1.transAxes,bbox=props)
ax1.legend()

ax1.set_xlabel('Time (UTC)')
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d %-H:%m"))
plt.setp( ax1.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax1.set_xlim(xlim)

ax2 = axs[2]

LHS_interp = np.interp(dt_list2[velmask],np.array([toTimestamp(time) for time in LHS_dt]),LHS)
ax2.plot(LHS_interp,RHS,linestyle='none',marker='o')
ax2.set_xlabel('LHS')
ax2.set_ylabel('RHS')
ax2.grid()
ax2.set_aspect(1)

# show plot
plt.subplots_adjust(hspace=0.5)
plt.show()