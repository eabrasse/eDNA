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

# drop NAN's from master list
dfm = dfm.dropna(subset=['lat', 'lon'])

# dolphin pen location
lon0 = -122.729779; lat0 = 47.742773

#DNA quants to master dataframe
# dfm['mean_est']=np.exp(dfq[dfq.data_ID.isin(dfm.data_ID)].mean_est.values)
# dfm['p0.25']=dfm['mean_est']-np.exp(dfq[dfq.data_ID.isin(dfm.data_ID)]['p0.25'].values)
# dfm['p0.75']=np.exp(dfq[dfq.data_ID.isin(dfm.data_ID)]['p0.75'].values)-dfm['mean_est']
dfm['mean_est']=dfq[dfq.data_ID.isin(dfm.data_ID)].mean_est.values
dfm['p0.25']=dfm['mean_est']-dfq[dfq.data_ID.isin(dfm.data_ID)]['p0.25'].values
dfm['p0.75']=dfq[dfq.data_ID.isin(dfm.data_ID)]['p0.75'].values-dfm['mean_est']
dfm['dist_from_pen'] = dfm.apply(lambda row: efun.ll2dist(row.lon,row.lat,lon0,lat0),axis=1) 

fig = plt.figure(figsize=(6, 4))

ax = fig.gca()

# ax.errorbar(dfm.dist_from_pen,dfm.mean_est,yerr=[dfm['p0.25'],dfm['p0.75']],linestyle='none',marker='o',mfc='None',mec='k',markersize=8)
ax.errorbar(dfm.dist_from_pen,dfm.mean_est,yerr=1.039,linestyle='none',marker='o',mfc='None',mec='k',markersize=8)

ax.set_xlabel('dist (m)')
ax.set_ylabel(r'log estimated DNA (copies $\mu \mathrm{L}^{-1}$)')
# ax.set_yscale('log')

x0 = dfm.dist_from_pen.values
y0 = dfm.mean_est.values
xinds = np.argsort(x0)
p = np.polyfit(x0[xinds],y0[xinds],1)
xx = np.linspace(x0.min(),x0.max(),100)
yy = xx*p[0]+p[1]

ax.plot(xx,yy,linestyle='dashed',color='m')
ax.text(0.5,0.7,'best fit line in log space = {:0.4f}x+{:0.4f}\n equiv in lin space = {:0.4f} exp({:0.4f}x)'.format(p[0],p[1],np.exp(p[1]),p[0]),transform=ax.transAxes,ha='center',va='center',color='m')

# fig.subplots_adjust(right=0.98,left=0.1,bottom=0.1,top = 0.98,wspace=0.5)
# outfn = home+'plots/Feb2023_moorings_obs.png'
# plt.savefig(outfn)
# print(f'Saved to {outfn}')
plt.show(block=False)
plt.pause(0.1)
# plt.close()


