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

# tide data: tide height (m) at tide gauge
tide_height = {}
tide_height['fn'] = home+'9445133_MLLW_20230131_20230202.txt'
tide_height['df'] = pd.read_csv(tide_height['fn'],sep='\s+',engine='python',skiprows=13)
tide_height['df']['datetime'] =  pd.to_datetime(tide_height['df']['Date']+' '+tide_height['df']['Time'])
tide_height['df']['datetime'] = tide_height['df']['datetime'].dt.tz_localize(tz='UTC')
tide_height['label'] = 'Tide height (m)'
# tide_height['df']['data'] = tide_height['df']['Pred']
# tide_height['df'].index=tide_height['df']['datetime']

tide_height_frac_change = {}
tide_height_frac_change['df'] = tide_height['df'].copy(deep=True)
tide_height_frac_change['label'] = 'Tide height fractional change (unitless)'
# tide_height_frac_change['df']['data'] = np.nan*tide_height['df']['data'] # "nans like tide_height['df']['data']"
tideheightchange = (tide_height['df']['Pred'][1:].values-tide_height['df']['Pred'][:-1].values)/(10+tide_height['df']['Pred'][:-1].values)
tideheightchange = 0.5*(tideheightchange[1:]+tideheightchange[:-1])
tideheightchange = np.pad(tideheightchange,(1,1),constant_values=np.nan)
tide_height_frac_change['df']['tideheightfracchange'] = tideheightchange
tide_height_frac_change['df'].index=tide_height_frac_change['df']['datetime']

over_height={}
over_height['df'] = tide_height['df'].copy(deep=True)
over_height['label'] = 'Over height (1/m)'
over_height['df']['overheight'] = 1/(10+tide_height['df']['Pred'].values)
over_height['df'].index=over_height['df']['datetime']



# # create dDNA/dt data
# dDNAdt = {}
# dDNAdt['df'] = pd.DataFrame()
# dDNAdt['label'] = 'dDNAdt (copies/uL/s)'
# myDNA = 'process'
# if myDNA=='ESP':
#     pt = ESP['df']['data0'].values
#     dDNA = np.diff(pt)
#     dt = np.diff(ESP['df']['datetime'].view('int64').values//10**9)
#     dDNAdt['df']['data0'] = dDNA/dt
#     pdt = ESP['df']['datetime'].values
#     dDNAdt['df']['datetime'] = [pdt[i]+np.timedelta64(int(dt[i]),'s') for i in range(len(pdt)-1)]
# else:
#     pt = process['df']['data0'].values
#     dDNA = np.diff(pt)
#     dt = np.diff(process['df']['datetime'].view('int64').values//10**9)
#     dDNAdt['df']['data0'] = dDNA/dt
#     pdt = process['df']['datetime'].values
#     dDNAdt['df']['datetime'] = [pdt[i]+np.timedelta64(int(dt[i]),'s') for i in range(len(pdt)-1)]

#~~~~~~~~~~~~~~~~~~~~
#
dataset_list = dolphin,tide_height_frac_change,over_height

# generate master list bookends while converting to timestamps
start_datetime_list = []
end_datetime_list = []
for dataset in dataset_list:
    df = dataset['df']
    # start_datetime_list.append(df.datetime.dt.floor('H').values[0])
    # end_datetime_list.append(df.datetime.dt.ceil('H').values[-1])
    df['timestamp'] = df.index.values.view('int64') // 10 ** 9 # convert to time stamp for easier comparisons and interpolation
    df.index=df['timestamp']

master = dolphin['df'].join(tide_height_frac_change['df'],lsuffix='_dolphin',rsuffix='_tideheightfracchange')
master = master.join(over_height['df'],rsuffix='_overheight')
master = master.dropna()
master = master[['ndolphin_net','tideheightfracchange','overheight']]


ESP_master = ESP['df'][['copies_per_mLseawater','timestamp_0','timestamp_1']]

# #~~~~~~~~~~~~~~~~~~~~~

# interpolate data onto master timestamps

plotting=False
saving=True
# if plotting:
#     fig = plt.figure(figsize=(8,8))
#     ax = fig.gca()
#     fig2 = plt.figure(figsize=(8,8))
#     ax2=fig2.gca()
for dataset in dataset_list:
    # note: not saving in dataframes bc it will have different time indexes!
    # dataset['data1'] = np.interp(np.array(master1['Timestamp (sec)']),np.array(dataset['df']['timestamp']),dataset['df']['data0'])
    # dataset['data'] = np.interp(np.array(master['Timestamp (sec)']),np.array(dataset['df']['timestamp']),dataset['df']['data0'])
    
    if plotting:
        fig = plt.figure(figsize=(8,8))
        ax = fig.gca()
        ax.plot(dataset['df'].index,dataset['df']['data'],linestyle='none',lw=1,color=tab20b(0),marker='o',markerfacecolor=tab20b(0),markeredgecolor=tab20b(0),label='raw data')
        # ax.plot(master['Datetime'],dataset['data'],linestyle='solid',lw=0.5,color=tab20b(5),marker='o',markerfacecolor='None',markeredgecolor=tab20b(5),label='interp to ESP time')
        # ax.plot(master1['Datetime'],dataset['data1'],linestyle='dashed',lw=0.5,color=tab20b(9),marker='x',markerfacecolor=tab20b(9),markeredgecolor=tab20b(9),label='interp to hourly time')
        ax.set_xlabel('Time')
        ax.set_ylabel('variable')
        ax.legend(loc='upper right')
        ax.text(0.1,0.9,dataset['label'],transform=ax.transAxes,ha='left',va='top')
        plt.show(block=False)
        # plt.pause(0.1)
fw,fh=efun.gen_plot_props()
fig = plt.figure(figsize=(fw*2,fh*2))
axlog = fig.add_subplot(2,1,1)
axlin = fig.add_subplot(2,1,2)
DNAcolor=tab20c(8)
dolphincolor='m'

ts0 = ESP['df']['timestamp_0'].min()
xscale = 1/3600
for ax in axlog,axlin:
    ax2 = ax.twinx()
    ax2.fill_between(xscale*(dolphin['df']['timestamp']-ts0),dolphin['df']['ndolphin_net'],color='none',edgecolor=dolphincolor,alpha=0.3,zorder=100,hatch='//')
    for i in ESP['df'].index:
        hatchstyle=''
        fillcolor= DNAcolor
        alp = 1.0
        # if ESP['df']['corr_conc'][i]>1000:
            # fillcolor=tab20c(8)
            # alp = 0.9
            # hatchstyle='\\'
        # else:
            # fillcolor=tab20c(11)
            # alp=0.4
            # hatchstyle=''
        ax.fill_between([xscale*(ESP['df']['timestamp_0'][i]-ts0),xscale*(ESP['df']['timestamp_1'][i]-ts0)],[ESP['df']['copies_per_mLseawater'][i]]*2,
            alpha=alp,color=fillcolor,zorder=5,edgecolor='none',hatch=hatchstyle)#,lw=4)

    
    # ax.set_yscale('log')
    # ax.set_ylabel(r'DNA concentration (copies $\mu \mathrm{L}^{-1}$)',color=DNAcolor)
    ax.set_ylabel(r'DNA concentration (copies $\mathrm{mL}^{-1}$)',color=DNAcolor)
    ax2.set_yticks([0,1,2,3])
    ax2.set_ylabel('N dolphins',color=dolphincolor)
    ax2.tick_params(axis='y',labelcolor=dolphincolor)
    ax.tick_params(axis='y',labelcolor=DNAcolor)
    fudge = 5000
    ax.set_xlim([xscale*(ESP['df']['timestamp_0'].min()-fudge-ts0),xscale*(ESP['df']['timestamp_1'].max()+fudge-ts0)])
    yl2 = ax2.get_ylim()
    ax2.set_ylim([0,yl2[1]])
    ax.axhline(y=1000,color=DNAcolor,linestyle='dashed',lw=1)

yl = axlin.get_ylim()
axlin.set_ylim([0,yl[1]])
axlog.text(1,1000,'Spike classification threshold',ha='left',va='bottom',color=DNAcolor,fontsize=8)
axlog.set_yscale('log')
axlin.set_xlabel('Hours since start of ESP sampling')

axlog.text(0.1,0.9,'a) Logarithmic scale',transform=axlog.transAxes)
axlin.text(0.1,0.9,'b) Linear scale',transform=axlin.transAxes)

fig.subplots_adjust(top=0.99,bottom=0.1)
plt.show(block=False)
plt.pause(0.1)
# outfn = '/Users/elizabethbrasseale/Projects/eDNA/plots/ESP_v_dolphin_ddPCR.png'
# plt.savefig(outfn)

# if plotting:
#     # ax.scatter(ESP['df']['datetime'],ESP['df']['PB_quantity'])
#
#     ax.set_yscale('log')
#     # plt.show(block=False)
#     # plt.pause(0.1)
#
#     ax2.plot(tide_height['df']['datetime'],tide_height['df']['data0'],color='blue',label='tide height')
#     ax2.plot(tide_height_frac_change['df']['datetime'],tide_height_frac_change['df']['data0'],color='magenta',label='fractional change')
#     ax2.set_xlabel('datetime')
#     ax2.legend()
#     ax2.grid()
#     plt.show(block=False)
#     plt.pause(0.1)

if saving:
    # master['DNA (copies/uL) (ESP triplicate mean)'] = ESP['data']
    # master['N (dolphins)'] = dolphin['data']
    # # master['Current velocity at tide gauge (m/s)'] = tide_gauge['data']
    # # master['Current velocity at dolphin pen (m/s)'] = tide_dolphinpen['data']
    # master['dDNAdt (copies/uL/s)'] = dDNAdt['data']
    # master['Tide height (m)'] = tide_height['data']
    # master['smoothed DNA (copies/uL) (observation error removed)'] = process['data']
    # master['Tide height fractional change (unitless)'] = tide_height_frac_change['data']
    # master['Over height (1/m)'] = over_height['data']
    # for dataset in dataset_list:
    #     master[dataset['label']] = dataset['data']


    outfn = home+'dolphin_tide_minutefreq.csv'
    master.to_csv(outfn)
    print('saved to {}'.format(outfn))
    
    outfn = home+'ESP_timestamps_mLseawater.csv'
    ESP_master.to_csv(outfn)
    print('saved to {}'.format(outfn))