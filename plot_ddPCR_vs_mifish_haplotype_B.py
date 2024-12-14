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
from scipy.stats import linregress

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
# 
ESP0={}
ESP0['fn'] = home+'MURI_Module1_datasheets - 03_ESP1_Feb2023.csv'
ESP0['df'] = pd.read_csv(ESP0['fn'],sep=',',engine='python')
ESP0['df'] = ESP0['df'].dropna(subset=['ESP_date','PB_quantity']) # ignore data without timestamps
ESP0['df']['datetime_0'] = pd.to_datetime(ESP0['df']['ESP_date']+' '+ESP0['df']['ESP_filter_t0'])
# ESP0['df']['datetime_1'] = pd.to_datetime(ESP0['df']['ESP_date']+' '+ESP0['df']['ESP_filter_t1'])
ESP0['df']['datetime_1'] = ESP0['df'].apply(lambda row: row.datetime_0+timedelta(hours=int(row['ESP_filter_totaltime'][0]),minutes=int(row['ESP_filter_totaltime'][2:4]),
    seconds=int(row['ESP_filter_totaltime'][5:])),axis=1)

ESP0['label'] = 'ESP raw data'

#then stamp those on triplicate means
ESP={}
# ESP['fn'] = home+'03_ESP1_Feb2023_hourly.csv'
ESP['fn'] = home+'ddPCR_data_2024-03-12_dPB105.csv'
ESP['df'] = pd.read_csv(ESP['fn'],sep=',',engine='python')
# ESP['df'] = ESP['df'].dropna()'
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
# ESP['df'].index=ESP['df'].datetime

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
#
dataset_list = [dolphin]


ESP_master = ESP['df'][['corr_conc','timestamp_0','timestamp_1']]

timestamp0_list = [t for t in ESP['df']['timestamp_0']]
timestamp1_list = [t for t in ESP['df']['timestamp_1']]

###
#explore mifish data

mifish = {}
mifish['fn'] = home+'ESP_ASVs_filtered.csv'
mifish['df'] = pd.read_csv(mifish['fn'],sep=',',engine='python')
mifish['df']['Sample code2'] = mifish['df']['Sample code'].str.slice(0,-5)
sample_list = [f'MFU.52345.{count:}' for count in range(171,213)]
# mifish
mifish['df']=mifish['df'][mifish['df']['Sample code2'].isin(sample_list)]
#combine types for first pass analysis
# mifish['df']=mifish['df0'].T.groupby(['_'.join(s.split('_')[:-1]) for s in mifish['df0'].T.index.values]).sum().T
mifish['df']['timestamp_0'] = timestamp0_list
mifish['df']['timestamp_1'] = timestamp1_list
column_names = list(mifish['df'].columns)
column_sub = [column for column in column_names if (column[:8]=='Tursiops')&(column[-2:]!='14')]
ndolphhap = len(column_sub)
total_Tursiops = np.zeros(len(mifish['df'].index))
for column in column_sub:
    total_Tursiops+=mifish['df'][column]

# mifish['df']['Tursiops'] = metamifish['df'][metamifish['df'].isin(column_sub)].sum()

# for column in column_sub:
#     mifish['df'][column] = mifish['df'][column].astype(float)
    
## add total reads to mifish
metamifish = {}
metamifish['fn'] = home+'ASV_table.csv'
metamifish['df'] = pd.read_csv(metamifish['fn'],sep=',',engine='python')
metamifish['df'] = metamifish['df'].groupby(metamifish['df'].Sample).sum()
mfmsample_list = [f'MFU-52345-{count:}' for count in range(171,213)]
metamifish['df'] = metamifish['df'][metamifish['df'].index.isin(mfmsample_list)]

running_min = 10
mifish['df']['nReads'] = metamifish['df']['nReads'].values
for column in column_sub:
    mifish['df'][column] = mifish['df'][column].astype(float)
    mifish['df'][column] = mifish['df'].apply(lambda row: row[column]/row.nReads,axis=1)
    mifish['df'][column][mifish['df'][column]==0] = 1*np.exp(-15)
###

fw,fh=efun.gen_plot_props()
fig = plt.figure(figsize=(fw*2,fh))
ax = fig.gca()
# DNAcolor=tab20c(8)
# DNAcolor=tab20c(11)
DNAcolor='k'
# dolphincolor='m'
dolphincolor=tab10(0)

ts0 = ESP['df']['timestamp_0'].min()
xscale = 1/3600


ax2 = ax.twinx()
ax.set_zorder(3)
ax2.set_zorder(1)
ax.set_frame_on(False)
# mifishcol=tab10(2)
mifish_cmap = plt.get_cmap('Set2',8)
ax2.plot(xscale*(0.5*(ESP['df']['timestamp_0']+ESP['df']['timestamp_1'])-ts0),np.log(np.array(ESP['df']['corr_conc'])),linestyle='--',marker='none',mfc='None',mec='None',color='k',lw=1,markersize=8)
# ax2.fill_between(xscale*(dolphin['df']['timestamp']-ts0),dolphin['df']['ndolphin_net'],color='none',edgecolor=dolphincolor,alpha=1,zorder=100,hatch='//')
# ax2.fill_between(xscale*(dolphin['df']['timestamp']-ts0),dolphin['df']['ndolphin_net'],color=dolphincolor,edgecolor='none',alpha=.5,zorder=1,)

hap_count=0
for column in column_sub:
    ax.plot(xscale*(0.5*(mifish['df']['timestamp_0']+mifish['df']['timestamp_1'])-ts0),np.log(np.array(mifish['df'][column])),linestyle='-',marker='o',mfc=mifish_cmap(hap_count),
    mec='None',color=mifish_cmap(hap_count),lw=1,markersize=8)
    
    hap_count+=1
    

# for i in mifish['df'].index:
#
#     running_sum = 0
#     hap_count=0
#     for column in column_sub:
#         # barheight=0
#         newcount = mifish['df'][column][i]
#         # if (newcount > 0):
#             # frac = newcount/total_Tursiops[i]
#             # barheight=frac*np.log(total_Tursiops[i])
#             # newcount = np.log(newcount)
#         ax.fill_between([xscale*(mifish['df']['timestamp_0'][i]-ts0),xscale*(mifish['df']['timestamp_1'][i]-ts0)],y1=[running_sum]*2,y2=[running_sum+newcount]*2,
#             alpha=1,color=mifish_cmap(hap_count),zorder=250,edgecolor='none')#,lw=4)
#         running_sum+=newcount
#         hap_count+=1

hap_count=0
for column in column_sub:
    ax.text(0.025,0.95-0.05*hap_count,column,color=mifish_cmap(hap_count),transform=ax.transAxes,ha='left',va='top',fontsize=10,fontweight='bold')
    hap_count+=1

ax2.set_ylabel(r'log DNA concentration (copies $\mathrm{mL}^{-1}$)',color='k')
# ax2.set_yticks([0,1,2,3])
ax.set_ylabel('log reads/Nreads',color=mifish_cmap(0))
ax2.tick_params(axis='y',labelcolor='k')
ax.tick_params(axis='y',labelcolor=mifish_cmap(0))
fudge = 5000
ax.set_xlim([xscale*(mifish['df']['timestamp_0'].min()-fudge-ts0),xscale*(mifish['df']['timestamp_1'].max()+fudge-ts0)])
# yl2 = ax2.get_ylim()
# ax2.set_ylim([0,yl2[1]])
ax.set_ylim([-9.8,-0.5])
# ax.axhline(y=20,color=DNAcolor,linestyle='dashed',lw=1)

# yl = axlin.get_ylim()
# axlin.set_ylim([0,yl[1]])
# ax.text(1,20,r'20 copies $\mathrm{mL}^{-1}$',ha='left',va='bottom',color=DNAcolor,fontsize=8)
# ax.set_yscale('log')
# ax2.set_yscale('log')
ax.set_xlabel('Hours since start of ESP sampling')

# axlog.text(0.1,0.9,'a) Logarithmic scale',transform=axlog.transAxes)
# axlin.text(0.1,0.9,'b) Linear scale',transform=axlin.transAxes)

fig.subplots_adjust(top=0.99,bottom=0.14)

# mifish['df']['Tursiops_truncatus'] = mifish['df'].Tursiops_truncatus.astype(float)
# fig2 = plt.figure(figsize=(fw*2,fh))
# ax2 = fig2.gca()
# logESP = np.log(np.array(ESP['df']['copies_per_mLseawater']))
# logmifish = np.log(np.array(mifish['df']['Tursiops_truncatus']))
# ax2.scatter(logESP,logmifish,c='k',s=20)
# ESP_inds = np.argsort(logESP)
# xx = logESP[ESP_inds]
# mifish_list = logmifish[ESP_inds]
# slope, intercept, r_value, p_value, std_err = linregress(xx,mifish_list)
# r_squared = r_value**2
# yy = intercept+slope*np.array(xx)
# ax2.plot(xx,yy,color=mifishcol,linestyle='dashed')
# ax2.text(0.9,0.1,'$\mathrm{R}^{2}$'+f' = {r_squared:.2f}',transform=ax2.transAxes,ha='right',color=mifishcol)
# ax2.set_xlabel('log ddPCR Tursiops DNA copies mL-1')
# ax2.set_ylabel('log # mifish Tursiops reads')

fig.show()

# fig2.subplots_adjust(top=0.99,bottom=0.14)
# fig2.show()
plt.pause(0.1)
# outfn = '/Users/elizabethbrasseale/Projects/eDNA/plots/ESP_ddPCR_vs_mifish.png'
# plt.savefig(outfn,format='png',dpi=400)
