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

timestamp_list = [t for t in ESP['df']['timestamp_0']]

###
#explore mifish data

mifish = {}
mifish['fn'] = home+'ESP_ASVs_filtered.csv'
mifish['df0'] = pd.read_csv(mifish['fn'],sep=',',engine='python')
mifish['df0']['Sample code2'] = mifish['df0']['Sample code'].str.slice(0,-5)
mfsample_list = [f'MFU.52345.{count:}' for count in range(171,213)]
# mifish
mifish['df0']=mifish['df0'][mifish['df0']['Sample code2'].isin(mfsample_list)]
#combine types for first pass analysis
mifish['df']=mifish['df0'].T.groupby(['_'.join(s.split('_')[:-1]) for s in mifish['df0'].T.index.values]).sum().T
mifish['df']['timestamp'] = timestamp_list
###

## add total reads to mifish
metamifish = {}
metamifish['fn'] = home+'ASV_table.csv'
metamifish['df'] = pd.read_csv(metamifish['fn'],sep=',',engine='python')
metamifish['df'] = metamifish['df'].groupby(metamifish['df'].Sample).sum()
mfmsample_list = [f'MFU-52345-{count:}' for count in range(171,213)]
metamifish['df'] = metamifish['df'][metamifish['df'].index.isin(mfmsample_list)]

mifish['df']['nReads'] = metamifish['df']['nReads'].values
column_names = list(mifish['df'].columns)
column_sub = [column for column in column_names if (column[:8]=='Tursiops')&(column[-2:]!='14')]
ndolphhap = len(column_sub)
for column in column_sub:
    mifish['df'][column] = mifish['df'][column].astype(float)
    mifish['df'][column] = mifish['df'].apply(lambda row: row[column]/row.nReads,axis=1)


fw,fh=efun.gen_plot_props()
fig = plt.figure(figsize=(fw*3,fh))
gs = GridSpec(1,3)
ax = fig.add_subplot(gs[:2])
# DNAcolor=tab20c(8)
# DNAcolor=tab20c(11)
DNAcolor='k'
# dolphincolor='m'
dolphincolor=tab10(0)

ts0 = ESP['df']['timestamp_0'].min()
xscale = 1/3600


ax2 = ax.twinx()
# ax.set_zorder(1)
# ax2.set_zorder(3)
# ax.set_frame_on(False)
mifishcol=tab10(2)
ax2.plot(xscale*(mifish['df']['timestamp']-ts0),mifish['df']['Tursiops_truncatus'],linestyle='-',marker='o',mfc=mifishcol,mec='None',color=mifishcol)
# ax2.fill_between(xscale*(dolphin['df']['timestamp']-ts0),dolphin['df']['ndolphin_net'],color='none',edgecolor=dolphincolor,alpha=1,zorder=100,hatch='//')
# ax2.fill_between(xscale*(dolphin['df']['timestamp']-ts0),dolphin['df']['ndolphin_net'],color=dolphincolor,edgecolor='none',alpha=.5,zorder=1,)
for i in ESP['df'].index:
    hatchstyle=''
    # if ESP['df']['corr_conc'][i]>1000:
    #     fillcolor=tab20c(8)
    # #     # alp = 0.9
    # #     # hatchstyle='\\'
    # else:
    #     fillcolor=tab20c(11)
    #     # alp=0.4
    #     # hatchstyle=''
    fillcolor= DNAcolor
    alp = 1.0
    ax.fill_between([xscale*(ESP['df']['timestamp_0'][i]-ts0),xscale*(ESP['df']['timestamp_1'][i]-ts0)],[ESP['df']['copies_per_mLseawater'][i]]*2,
        alpha=alp,color=fillcolor,zorder=250,edgecolor='none',hatch=hatchstyle)#,lw=4)
    
    # if i>2:
        # ax.plot([xscale*(ESP['df']['timestamp_0'][i-1]-ts0),xscale*(ESP['df']['timestamp_1'][i-1]-ts0),xscale*(ESP['df']['timestamp_0'][i]-ts0),xscale*(ESP['df']['timestamp_1'][i]-ts0)],
        # [ESP['df']['corr_conc'][i-1],ESP['df']['corr_conc'][i-1],ESP['df']['corr_conc'][i],ESP['df']['corr_conc'][i]],color=DNAcolor,linewidth=0.5,zorder=500,linestyle='solid')



ax.set_ylabel(r'DNA concentration (copies $\mathrm{mL}^{-1}$)',color=DNAcolor)
# ax2.set_yticks([0,1,2,3])
ax2.set_ylabel('Reads/nReads',color=mifishcol)
ax2.tick_params(axis='y',labelcolor=mifishcol)
ax.tick_params(axis='y',labelcolor=DNAcolor)
fudge = 5000
ax.set_xlim([xscale*(ESP['df']['timestamp_0'].min()-fudge-ts0),xscale*(ESP['df']['timestamp_1'].max()+fudge-ts0)])
# yl2 = ax2.get_ylim()
# ax2.set_ylim([0,yl2[1]])
# ax.axhline(y=20,color=DNAcolor,linestyle='dashed',lw=1)

# yl = axlin.get_ylim()
# axlin.set_ylim([0,yl[1]])
# ax.text(1,20,r'20 copies $\mathrm{mL}^{-1}$',ha='left',va='bottom',color=DNAcolor,fontsize=8)
ax.set_yscale('log')
ax2.set_yscale('log')
ax.set_xlabel('Hours since start of ESP sampling')

ax.text(0.1,0.9,'a)',transform=ax.transAxes)
# axlog.text(0.1,0.9,'a) Logarithmic scale',transform=axlog.transAxes)
# axlin.text(0.1,0.9,'b) Linear scale',transform=axlin.transAxes)

# fig.subplots_adjust(top=0.99,bottom=0.14)

mifish['df']['Tursiops_truncatus'] = mifish['df'].Tursiops_truncatus.astype(float)
# fig2 = plt.figure(figsize=(fw*2,fh))
ax2 = fig.add_subplot(gs[2])
logESP = np.log(np.array(ESP['df']['copies_per_mLseawater']))
logmifish = np.log(np.array(mifish['df']['Tursiops_truncatus']))
ax2.scatter(logESP,logmifish,c='k',s=20)
ESP_inds = np.argsort(logESP)
xx = logESP[ESP_inds]
mifish_list = logmifish[ESP_inds]
slope, intercept, r_value, p_value, std_err = linregress(xx,mifish_list)
r_squared = r_value**2
yy = intercept+slope*np.array(xx)
ax2.plot(xx,yy,color=mifishcol,linestyle='dashed')
ax2.text(0.9,0.1,'$\mathrm{R}^{2}$'+f' = {r_squared:.2f}',transform=ax2.transAxes,ha='right',color=mifishcol)
ax2.set_xlabel('log DNA copies mL-1')
ax2.set_ylabel('log Reads/nReads')
yl = ax2.get_ylim()
xl = ax2.get_xlim()
aspect = (xl[1]-xl[0])/(yl[1]-yl[0])
ax2.set_aspect(aspect)

ax2.text(0.1,0.9,'b)',transform=ax2.transAxes)

fig.subplots_adjust(top=0.99,bottom=0.14,wspace=0.8,right=0.98,left=0.1)
# fig.show()
# fig2.show()
# plt.pause(0.1)
outfn = '/Users/elizabethbrasseale/Projects/eDNA/plots/ESP_ddPCR_vs_mifish_B.png'
plt.savefig(outfn,format='png',dpi=400)
