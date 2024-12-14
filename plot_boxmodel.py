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

plt.close('all')
tab10 = plt.get_cmap('tab10',10)
tz = pytz.utc
pdt = pytz.timezone('America/Vancouver')

# load in data
home = '/Users/elizabethbrasseale/Projects/eDNA/data/'


ESP_fn = home+'MURI_Module1_datasheets - 03_ESP1_Feb2023.csv'

ESP_df = pd.read_csv(ESP_fn,sep=',',engine='python')
ESP_df['datetime'] = pd.to_datetime(ESP_df['ESP_date']+' '+ESP_df['ESP_filter_t0'])
ESP_df = ESP_df.dropna(subset=['ESP_date','PB_quantity'])
ESP_df['datetime'] = ESP_df['datetime'].dt.tz_localize(tz='America/Vancouver')
ESP_df['datetime'] = ESP_df['datetime'].dt.tz_convert('UTC')
datetimes = ESP_df['datetime'].unique()

ddPCR_fn = home+'ddPCR_data_2024-03-12_dPB105.csv'
ddPCR_df = pd.read_csv(ddPCR_fn,sep=',',engine='python')
samples_to_drop = ['52345-169','52345-170','52345-213','52345-214','52345-215','52345-216','52345-217','52345-218','52345-219','NTC','offtarget']
ddPCR_df = ddPCR_df[~ddPCR_df['sample'].isin(samples_to_drop)]

datafn = home+'ESP_box_model_output.csv'
df = pd.read_csv(datafn,sep=',',engine='python')
df['Datetime']=pd.to_datetime(df['Datetime']).dt.tz_localize(tz='UTC')

fw,fh = efun.gen_plot_props()
fig = plt.figure(figsize=(7,5))
ax0A = fig.gca()

# ax.plot(ESP_df['datetime'],ESP_df['PB_quantity'],linestyle='None',marker='o',mec=tab10(0),mfc = 'None')
ax0A.plot(datetimes,ddPCR_df['corr_conc'],linestyle='None',marker='o',mec=tab10(0),mfc = tab10(0))
ax0B = ax0A.twinx()
ax0B.plot(df['Datetime'],df['N..dolphins.'],linestyle='solid',color=tab10(2))
ax0B.set_ylabel('N Dolphins',color=tab10(2))
ax0B.set_yticks([0,1,2,3])
# ax.text(0.02,0.98,r'Observations, $y$',color=tab10(0),transform=ax.transAxes,ha='left',va='top')
# ax.plot(df['Datetime'],df['modelDNAconc'],linestyle='dashed',marker='None',color=tab10(1))
# ax.text(0.02,0.93,r'Model, $y$~N($\mu$,$\sigma^{2}$)'+'\n'+r'     $\mu_{t}=\mu_{t-1}(\alpha - \frac{\Delta SSH}{H})+e^{\beta \lambda}$'+'\n'+r'     $\lambda$~Pois(Dolph)',transform=ax.transAxes,ha='left',va='top',color=tab10(1))
# slope, intercept, r_value, p_value, std_err = linregress(df['modelDNAconc'],df['DNA..copies.uL...ESP.triplicate.mean.'])

ax0A.set_xlabel('Date, Time (UTC)')
ax0A.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d, %H:%M"))
plt.setp( ax0A.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
# ax.set_yscale('log')
ax0A.set_ylabel(r'DNA conc (copies $\mu\mathrm{L}^{-1}$)',color=tab10(0))


fig2 = plt.figure(figsize=(7,5))
ax1A = fig2.gca()
ax1A.plot(datetimes,ddPCR_df['corr_conc'],linestyle='None',marker='o',mec=tab10(0),mfc = tab10(0))
ax1B = ax1A.twinx()
ax1B.plot(df['Datetime'],df['Tide.height..m.'],linestyle='solid',color=tab10(4))
ax1B.set_ylabel('SSH (m)',color=tab10(4))
# ax3.set_yticks([0,1,2,3])

ax1A.set_xlabel('Date, Time (UTC)')
ax1A.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d, %H:%M"))
plt.setp( ax1A.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
# ax.set_yscale('log')
ax1A.set_ylabel(r'DNA conc (copies $\mu\mathrm{L}^{-1}$)',color=tab10(0))

for ax in ax0A,ax1A:
    ax.set_yscale('log')
fig.subplots_adjust(left=0.15,bottom=0.23,right=0.87,top=0.99)
fig2.subplots_adjust(left=0.15,bottom=0.23,right=0.87,top=0.99)


# plt.subplots_adjust(left=0.15,bottom=0.23,right=0.87,top=0.99)
plt.show(block=False)
plt.pause(0.1)
