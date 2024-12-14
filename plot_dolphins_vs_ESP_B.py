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

plt.close('all')

tz = pytz.utc
pdt = pytz.timezone('America/Vancouver')
tab20c = plt.get_cmap('tab20c',20)
tab10 = plt.get_cmap('tab10',10)

# load in data
home = '/Users/elizabethbrasseale/Projects/eDNA/data/'

# ESP data with observation error removed
vol_extract = 100 # muL

# To get timestamps for process term, go back to raw ESP data
ESP0={}
ESP0['fn'] = home+'MURI_Module1_datasheets - 03_ESP1_Feb2023.csv'
ESP0['df'] = pd.read_csv(ESP0['fn'],sep=',',engine='python')
ESP0['df'] = ESP0['df'].dropna(subset=['ESP_date','PB_quantity']) # ignore data without timestamps
ESP0['df']['datetime_0'] = pd.to_datetime(ESP0['df']['ESP_date']+' '+ESP0['df']['ESP_filter_t0'])
ESP0['df']['datetime_1'] = ESP0['df'].apply(lambda row: row.datetime_0+timedelta(hours=int(row['ESP_filter_totaltime'][0]),minutes=int(row['ESP_filter_totaltime'][2:4]),
    seconds=int(row['ESP_filter_totaltime'][5:])),axis=1)

ESP0['label'] = 'ESP raw data'

#then stamp those on triplicate means
ESP={}
ESP['fn'] = home+'ddPCR_data_2024-03-12_dPB105.csv'
ESP['df'] = pd.read_csv(ESP['fn'],sep=',',engine='python')
sample_list = [f'52345-{count:}' for count in range(171,213)]
ESP['df']=ESP['df'][ESP['df']['sample'].isin(sample_list)]
ESP['df']['datetime_0']= ESP0['df']['datetime_0'].unique() # remove replicates of timestamps
ESP['df']['datetime_1']= ESP0['df']['datetime_1'].unique() # remove replicates of timestamps
ESP['df']=ESP['df'].sort_values(by='datetime_0')
ESP['df']['datetime_0'] = ESP['df']['datetime_0'].dt.tz_localize(tz='America/Vancouver')
ESP['df']['datetime_0'] = ESP['df']['datetime_0'].dt.tz_convert('UTC')
ESP['df']['timestamp_0'] = ESP['df']['datetime_0'].view('int64').values//10**9
ESP['df']['datetime_1'] = ESP['df']['datetime_1'].dt.tz_localize(tz='America/Vancouver')
ESP['df']['datetime_1'] = ESP['df']['datetime_1'].dt.tz_convert('UTC')
ESP['df']['timestamp_1'] = ESP['df']['datetime_1'].view('int64').values//10**9
ESP['label'] = 'DNA (copies/uL) (ESP triplicate mean)'
ESP['df']['copies_per_ext'] = ESP['df']['corr_conc']
ESP['df']['vol_filtered'] = ESP0['df']['filter_vol_ml'].copy()
ESP['df']['copies_per_mLseawater'] = ESP['df'].apply(lambda row: row.copies_per_ext*vol_extract/row.vol_filtered,axis=1)

# dolphin occupancy data, digitized from Navy clipboard
dolphin0 = {}
dolphin0['fn'] = home+'dolphin_presence.csv'
dolphin0['df'] = pd.read_csv(dolphin0['fn'],sep=',',engine='python')
dolphin0['df'] = dolphin0['df'].dropna()
dolphin0['df']['datetime'] = pd.to_datetime(dolphin0['df'].date+' '+dolphin0['df'].time)
dolphin0['df']['datetime'] = dolphin0['df']['datetime'].dt.tz_localize(tz='America/Vancouver')
dolphin0['df']['datetime'] = dolphin0['df']['datetime'].dt.tz_convert('UTC')
dolphin0['label'] = 'N (dolphins)'
dolphin0['df']['data'] = dolphin0['df']['ndolphin_net']
dolphin0['df'].index=dolphin0['df']['datetime']

dolphin = {}
dolphin['df'] = dolphin0['df'].resample('1min').ffill() # put dolphin data on 1m time series
dolphin['label'] = 'N (dolphins)'
dolphin['df']['timestamp'] = dolphin['df'].index.values.view('int64') // 10 ** 9 # convert to time stamp for easier comparisons and interpolation
dolphin['df'].index=dolphin['df']['timestamp']

#open figure
fw,fh=efun.gen_plot_props()
fig = plt.figure(figsize=(fw*2,fh))
ax = fig.gca()

# set some colors here
DNAcolor='k'
dolphincolor=tab10(0)
ts0 = ESP['df']['timestamp_0'].min() #doing this to transform x-axis relative to ESP start time
xscale = 1/3600 #then divide time stamp in seconds by 3600 into "hours since ESP start"

# since the dolphin occupancy data is on a very different magnitude scale,
# open a second y-axis that shares the same x-axis. The yticks for the dolphin data will be on the RHS
ax2 = ax.twinx()
ax.set_zorder(1)
ax2.set_zorder(3)
ax.set_frame_on(False) # some bullshit to make the transparent fill show up BENEATH the ESP bars instead of on top of it

# plot!
# dolphin data fill
ax2.fill_between(xscale*(dolphin['df']['timestamp']-ts0),dolphin['df']['ndolphin_net'],color=dolphincolor,edgecolor='none',alpha=.5,zorder=1)
# since 
for i in ESP['df'].index:
    ax.fill_between([xscale*(ESP['df']['timestamp_0'][i]-ts0),xscale*(ESP['df']['timestamp_1'][i]-ts0)],[ESP['df']['copies_per_mLseawater'][i]]*2,
        alpha=1.0,color=DNAcolor,zorder=250,edgecolor='none')
    

ax.set_ylabel(r'DNA concentration (copies $\mathrm{mL}^{-1}$)',color=DNAcolor)
ax2.set_yticks([0,1,2,3])
ax2.set_ylabel('N (dolphins)',color=dolphincolor)
ax2.tick_params(axis='y',labelcolor=dolphincolor)
ax.tick_params(axis='y',labelcolor=DNAcolor)
fudge = 5000
ax.set_xlim([xscale*(ESP['df']['timestamp_0'].min()-fudge-ts0),xscale*(ESP['df']['timestamp_1'].max()+fudge-ts0)])
yl2 = ax2.get_ylim()
ax2.set_ylim([0,yl2[1]])
ax.axhline(y=20,color=DNAcolor,linestyle='dashed',lw=1)


ax.text(1,20,r'20 copies $\mathrm{mL}^{-1}$',ha='left',va='bottom',color=DNAcolor,fontsize=8)
ax.set_yscale('log')
ax.set_xlabel('Hours since start of ESP sampling')

fig.subplots_adjust(top=0.99,bottom=0.14)
plt.show(block=False)
plt.pause(0.1)
# outfn = '/Users/elizabethbrasseale/Projects/eDNA/plots/ESP_v_dolphin_ddPCR_mLseawater.png'
# plt.savefig(outfn,format='png',dpi=400)
