import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from PIL import Image
import numpy as np
from datetime import datetime, timedelta
import pytz
import netCDF4 as nc
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.dates as mdates
import pickle
import string
import efun
import pytz
from scipy import stats

tz = pytz.utc
pdt = pytz.timezone('America/Vancouver')

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)
    

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'

atoz = string.ascii_lowercase

tab10 = plt.get_cmap('tab10',10)
Paired = plt.get_cmap('Paired',12)


#use Ryan's estimated eDNA quants
DNA_quant_fn = home+'data/projected_logDNA_concentrations.csv'
dfq = pd.read_csv(DNA_quant_fn,sep=',',engine='python')

# get metadata from master spreadsheet
master_fn = home+'data/All_Mod1_Molecular_Data - all_data.csv'
dfm = pd.read_csv(master_fn,sep=',',engine='python')

#create a unique identifier per sample including lab ID, biorep, and technrep
for df in dfq,dfm:
    df['data_ID'] = df.apply(lambda row: row.lab_ID + str(row.bio_rep) + str(row.tech_rep),axis=1)

# just look at matrix samples
dfm = dfm[dfm['task']=='matrix']
dfm = dfm[~dfm['reference'].str.startswith('SR_later_PCI_')]
# dfm.loc[dfm['collection_time']=='18:30','collection_time']='15:00' # fix inaccurate value in spreadsheet
dfm = dfm[dfm['collection_time']!='18:30']
#DNA quants to master dataframe
dfm['mean_est']=dfq[dfq.data_ID.isin(dfm.data_ID)].mean_est.values
dfm['p0.25']=dfm['mean_est']-dfq[dfq.data_ID.isin(dfm.data_ID)]['p0.25'].values
dfm['p0.75']=dfq[dfq.data_ID.isin(dfm.data_ID)]['p0.75'].values-dfm['mean_est']
# dfm['dist_from_pen'] = dfm.apply(lambda row: efun.ll2dist(row.lon,row.lat,lon0,lat0),axis=1) 
dfm['hour'] = dfm.apply(lambda row: row.collection_time[:2],axis=1)
dfm['minute'] = dfm.apply(lambda row: row.collection_time[3:],axis=1)
# dfm['datetime'] = pd.to_datetime('{}-{}-{} {}:00'.format(str(dfm.year),str(dfm.month),str(dfm.day),str(dfm.collection_time)))
dfm['datetime'] = pd.to_datetime({'year':dfm.year,'month':dfm.month,'day':dfm.day,'hour':dfm.hour,'minute':dfm.minute})
dfm['timestamp'] = dfm['datetime'].astype(int)
dfm['seconds_since_t0'] = dfm.timestamp-dfm.timestamp.values[0]
dfm['seconds_since_t0'] = dfm.seconds_since_t0/10**9

fig = plt.figure(figsize=(6, 4))

ax = fig.gca()

ax.plot(dfm.datetime,dfm.dolphin_DNA_copy_uL,linestyle='none',marker='o',mfc='None',mec=tab10(0),markersize=8)
# ax.plot(dfm.datetime,dfm.mean_est,linestyle='none',marker='o',mfc='None',mec='k',markersize=8)

# ax.errorbar(dfm.datetime,dfm.mean_est,yerr=[dfm['p0.25'],dfm['p0.75']],linestyle='none',marker='o',mfc='None',mec='k',markersize=8)
# ax.errorbar(dfm.datetime,dfm.mean_est,yerr=1.039,linestyle='none',marker='o',mfc='None',mec='k',markersize=8)

ax.set_xlabel('Time (PDT)')
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.set_ylabel(r'DNA (copies $\mu \mathrm{L}^{-1}$)')
ax.set_ylim([2000.05552794096448, 12643200.66871178932])
ax.set_yscale('log')

x0 = dfm.seconds_since_t0.values
y0 = np.log(dfm.dolphin_DNA_copy_uL.values)
xinds = np.argsort(x0)
# p = np.polyfit(x0[xinds],y0[xinds],1)
slope, intercept, r_value, p_value, std_err = stats.linregress(x0[xinds],y0[xinds])
xx = np.linspace(x0.min(),x0.max(),100)
yy = np.exp(intercept)*np.exp(xx*slope)
xx2plot = pd.date_range(dfm.datetime.min(),dfm.datetime.max(),100)

ax.plot(xx2plot,yy,linestyle='dashed',color=tab10(1))
# ax.text(0.5,0.6,'best fit line = {:0.2E}exp(-{:0.2E}x)'.format(np.exp(intercept),np.abs(slope)),transform=ax.transAxes,ha='center',va='center',color=tab10(1))

fig.subplots_adjust(right=0.98,left=0.2,bottom=0.15,top = 0.98)
# outfn = home+'plots/Feb2023_moorings_obs.png'
# plt.savefig(outfn)
# print(f'Saved to {outfn}')
plt.show(block=False)
plt.pause(0.1)
# plt.close()


