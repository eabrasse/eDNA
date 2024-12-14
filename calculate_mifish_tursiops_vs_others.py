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
import re

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
old_column_names = mifish['df0'].columns
new_column_names = [re.sub(r"\_"," ",re.sub(r"\d+","",col_name)).strip() for col_name in old_column_names]
mifish['df0'].columns=new_column_names
# mifish['df']=mifish['df0'].T.groupby(['_'.join(s.split('_')[:-1]) for s in mifish['df0'].T.index.values]).sum().T
mifish['df'] = mifish['df0'].T.groupby(new_column_names).sum().T
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
mifish['df']['Tursiops truncatus'] = mifish['df']['Tursiops truncatus'].astype(float)
column_names = list(mifish['df'].columns)
column_sub = [column for column in column_names if column not in ['Sample code','Sample code2','timestamp','nReads','','NA']]
ndolphhap = len(column_sub)
Tcorr = np.zeros((ndolphhap))
i=0
for column in column_sub:
    mifish['df'][column] = mifish['df'][column].astype(float)
    mifish['df'][column] = mifish['df'].apply(lambda row: row[column]/row.nReads,axis=1)
    
    Tcorr[i] = np.corrcoef(mifish['df']['Tursiops truncatus'],mifish['df'][column])[0,1]
    mifish['df'][column][mifish['df'][column]==0] = 10**-15
    i+=1
D = pd.DataFrame(data=Tcorr,index=column_sub)
D = D.sort_values(0)[-10:-1][::-1]
Ttcol = tab10(1)
Clpcol = tab10(4)
ts0 = ESP['df']['timestamp_0'].min()
xscale = 1/3600

fw,fh=efun.gen_plot_props()
fig = plt.figure(figsize=(fw*4,fh*2))
gs = GridSpec(3,3)

count=0
for spp in D.index:
    row= int(np.floor(count/3))
    col = int(np.mod(count,3))
    ax = fig.add_subplot(gs[row,col])
    # ax2 = ax.twinx()

    ax.plot(xscale*(mifish['df']['timestamp']-ts0),np.log(mifish['df']['Tursiops truncatus']),linestyle='-',marker='o',mfc=Ttcol,mec='None',color=Ttcol)
    inds = np.argwhere(np.log(mifish['df']['Tursiops truncatus'])>-2)
    for ind in inds:
        xi = ind[0]
        ax.fill_between(xscale*(mifish['df']['timestamp'][xi-1:xi+2]-ts0),[-15]*3,[0]*3,color=Ttcol,alpha=0.15)
    # ax.plot(xscale*(mifish['df']['timestamp'][inds]-ts0),np.log(mifish['df']['Tursiops truncatus'][inds]),linestyle='none',marker='*',markersize=12,mfc='yellow')
    ax.plot(xscale*(mifish['df']['timestamp']-ts0),np.log(mifish['df'][spp]),linestyle='-',marker='o',mfc=Clpcol,mec='None',color=Clpcol)
    if count==0:
        ax.text(0.1,0.9,'Tursiops truncatus',transform=ax.transAxes,color=Ttcol)
        ax.text(0.1,0.8,f'{spp} (R={D.values[count][0]:.2})', transform=ax.transAxes,color=Clpcol)
    else:
        ax.text(0.1,0.9,f'{spp} (R={D.values[count][0]:.2})', transform=ax.transAxes,color=Clpcol)
    
    ax.set_ylim([-10,0])

    # ax.set_yscale('log')
    # ax2.set_yscale('log')
    if col==0:
        ax.set_ylabel(r'log Reads/nReads')
        # ax.tick_params(axis='y',labelcolor=Ttcol)
    else:
        ax.set_yticklabels([''])
        
    if row<2:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('Hours since start of ESP sampling')
    # ax2.set_yticks([0,1,2,3])
    # ax2.set_ylabel(f'log {spp:} Reads/nReads',color=Clpcol)
    # ax2.tick_params(axis='y',labelcolor=Clpcol)
    
    fudge = 5000
    ax.set_xlim([xscale*(ESP['df']['timestamp_0'].min()-fudge-ts0),xscale*(ESP['df']['timestamp_1'].max()+fudge-ts0)])
    # yl2 = ax2.get_ylim()
    # ax2.set_ylim([0,yl2[1]])
    # ax.axhline(y=20,color=DNAcolor,linestyle='dashed',lw=1)

    # yl = axlin.get_ylim()
    # axlin.set_ylim([0,yl[1]])
    # ax.text(1,20,r'20 copies $\mathrm{mL}^{-1}$',ha='left',va='bottom',color=DNAcolor,fontsize=8)

    # 

    # axlog.text(0.1,0.9,'a) Logarithmic scale',transform=axlog.transAxes)
    # axlin.text(0.1,0.9,'b) Linear scale',transform=axlin.transAxes)
    count+=1

fig.subplots_adjust(top=0.99,bottom=0.1,right=0.9,left=0.1,hspace=0.1,wspace=0.1)


# fig.show()

# plt.pause(0.1)
outfn = '/Users/elizabethbrasseale/Projects/eDNA/plots/Tt_vs_top_R_ASVs.png'
plt.savefig(outfn,format='png',dpi=400)
