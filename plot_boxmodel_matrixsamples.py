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
import efun
# from scipy.stats import linregress
import matplotlib.dates as mdates
from matplotlib.dates import DayLocator, HourLocator, MonthLocator,DateFormatter, drange

plt.close('all')
tab10 = plt.get_cmap('tab10',10)
tz = pytz.utc
pdt = pytz.timezone('America/Vancouver')

# load in data
home = '/Users/elizabethbrasseale/Projects/eDNA/'


ESP_fn = home+'data/MURI_Module1_datasheets - 03_ESP1_Feb2023.csv'
ESP_df = pd.read_csv(ESP_fn,sep=',',engine='python')
ESP_df['datetime'] = pd.to_datetime(ESP_df['ESP_date']+' '+ESP_df['ESP_filter_t0'])
ESP_df = ESP_df.dropna(subset=['ESP_date','PB_quantity'])
ESP_df['datetime'] = ESP_df['datetime'].dt.tz_localize(tz='America/Vancouver')
ESP_df['datetime'] = ESP_df['datetime'].dt.tz_convert('UTC')

datafn = home+'data/ESP_box_model_output.csv'
df = pd.read_csv(datafn,sep=',',engine='python')
df['Datetime']=pd.to_datetime(df['Datetime']).dt.tz_localize(tz='UTC')

fw,fh = efun.gen_plot_props()
# fig, (ax2, ax) = plt.subplots(1,2, sharey=True,figsize=(fw*2,fh))
fig = plt.figure(figsize=(fw*2,fh))
# ax = fig.gca()
gs = GridSpec(1,3)
ax2 = fig.add_subplot(gs[0])
ax = fig.add_subplot(gs[1:], sharey=ax2)

ax.plot(ESP_df['datetime'],ESP_df['PB_quantity'],linestyle='None',marker='o',mec=tab10(0),mfc = 'None',markersize=8)
ax.text(0.02,0.98,r'Observations, $y$',color=tab10(0),transform=ax.transAxes,ha='left',va='top')
ax.plot(df['Datetime'],df['modelDNAconc'],linestyle='dashed',marker='None',color=tab10(1))
ax.text(0.02,0.93,r'Model, $y$~N($\mu$,$\sigma^{2}$)'+'\n'+r'     $\mu_{t}=\mu_{t-1}(\alpha - \frac{\Delta SSH}{H})+e^{\beta \lambda}$'+'\n'+r'     $\lambda$~Pois(Dolph)',transform=ax.transAxes,ha='left',va='top',color=tab10(1))
# slope, intercept, r_value, p_value, std_err = linregress(df['modelDNAconc'],df['DNA..copies.uL...ESP.triplicate.mean.'])

ax.set_xlabel('Date, Time (UTC)',ha='center', x=0.2)
for ax0 in ax,ax2:
    ax0.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d, %H:%M"))
    plt.setp( ax0.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax.set_yscale('log')


## add matrix sample plotting code

#use Ryan's estimated eDNA quants
DNA_quant_fn = home+'data/projected_logDNA_concentrations.csv'
dfq = pd.read_csv(DNA_quant_fn,sep=',',engine='python')

# get metadata from master spreadsheet
master_fn = home+'data/All_Mod1_Molecular_Data - all_data.csv'
dfm = pd.read_csv(master_fn,sep=',',engine='python')

#create a unique identifier per sample including lab ID, biorep, and technrep
for df0 in dfq,dfm:
    df0['data_ID'] = df0.apply(lambda row: row.lab_ID + str(row.bio_rep) + str(row.tech_rep),axis=1)

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

ax2.plot(dfm.datetime,dfm.dolphin_DNA_copy_uL,linestyle='none',marker='o',mfc='None',mec=tab10(0),markersize=8)
ax2.set_ylabel(r'DNA conc (copies $\mu\mathrm{L}^{-1}$)')
ax2.xaxis.set_major_locator(mdates.MinuteLocator(range(0,60,10)))

# do some formatting
ax2.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.tick_params(which='both',left=False)
ax.yaxis.tick_right()

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-0.5*d, 0.5*d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot( (-0.5*d, +0.5*d), (1 - d, 1 + d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((1-d, 1+d), (-d, +d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax2.text(0.05,0.05,'Inside heated pen',va='center',ha='left',color='gray',transform=ax2.transAxes)
ax.text(0.95,0.05,'Outside net pen',va='center',ha='right',color='gray',transform=ax.transAxes)

plt.subplots_adjust(left=0.1,bottom=0.23,right=0.9,top=0.95,wspace=0.05)
plt.show(block=False)
plt.pause(0.1)
# outfn = home+ 'plots/ESP_and_matrix_samples_B.png'
# plt.savefig(outfn,format='png',dpi=600)