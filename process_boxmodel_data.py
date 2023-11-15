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

plt.close('all')

tz = pytz.utc
pdt = pytz.timezone('America/Vancouver')

# load in data
home = '/Users/elizabethbrasseale/Projects/eDNA/data/'

# ESP data with observation error removed
# Doesn't have timestamps!
process = {}
process['fn'] = home+'process_term_mean.csv'
process['df'] = pd.read_csv(process['fn'],sep=',',engine='python')

# To get timestamps for process term, go back to raw ESP data
ESP = {}
ESP['fn'] = home+'MURI_Module1_datasheets - 03_ESP1_Feb2023.csv'
ESP['df'] = pd.read_csv(ESP['fn'],sep=',',engine='python')
ESP['df']['datetime'] = pd.to_datetime(ESP['df']['ESP_date']+' '+ESP['df']['ESP_filter_t0'])
ESP['df'] = ESP['df'].dropna(subset=['ESP_date','PB_quantity']) # ignore data without timestamps
process['df']['datetime']= ESP['df']['datetime'].unique() # remove replicates of timestamps
process['df']=process['df'].sort_values(by='datetime')
process['df']['datetime'] = process['df']['datetime'].dt.tz_localize(tz='America/Vancouver')
process['df']['datetime'] = process['df']['datetime'].dt.tz_convert('UTC')

# include ESP triplicate averaged
ESP_trip_avg = {}
ESP_trip_avg['fn'] = home+'03_ESP1_Feb2023_hourly.csv'
ESP_trip_avg['df'] = pd.read_csv(ESP_trip_avg['fn'],sep=',',engine='python')
ESP_trip_avg['df'] = ESP_trip_avg['df'].dropna()
ESP_trip_avg['df']['datetime'] = process['df']['datetime'].copy()

# dolphin occupancy data, digitized from Navy clipboard
dolphin = {}
dolphin['fn'] = home+'dolphin_presence.csv'
dolphin['df'] = pd.read_csv(dolphin['fn'],sep=',',engine='python')
dolphin['df'] = dolphin['df'].dropna()
dolphin['df']['datetime'] = pd.to_datetime(dolphin['df'].date+' '+dolphin['df'].time)
dolphin['df']['datetime'] = dolphin['df']['datetime'].dt.tz_localize(tz='America/Vancouver')
dolphin['df']['datetime'] = dolphin['df']['datetime'].dt.tz_convert('UTC')

tidecurrent = 'velocity'

if tidecurrent=='speed':
    # tide data: tidal speed = sqrt(u^2+v^2) extracted from model at tide gauge point and time stamp in UTC
    tide_gauge = {}
    tide_gauge['fn'] = home+'HC_dolph_tidal_current_ATTideGauge_9445133.csv'
    tide_gauge['df'] = pd.read_csv(tide_gauge['fn'],sep=',',engine='python')
    tide_gauge['df']['datetime'] = pd.to_datetime(tide_gauge['df']['Date Time'])
    tide_gauge['df']['datetime'] = tide_gauge['df']['datetime'].dt.tz_localize(tz='UTC')

    # tide data: tidal speed = sqrt(u^2+v^2) extracted from model at tide gauge point and time stamp in UTC
    tide_dolphinpen = {}
    tide_dolphinpen['fn'] = home+'HC_dolph_tidal_current_ATBangorDolphinPen.csv'
    tide_dolphinpen['df'] = pd.read_csv(tide_dolphinpen['fn'],sep=',',engine='python')
    tide_dolphinpen['df']['datetime'] = pd.to_datetime(tide_dolphinpen['df']['datetime'])
    tide_dolphinpen['df']['datetime'] = tide_dolphinpen['df']['datetime'].dt.tz_localize(tz='UTC')
    
else:
    # tide data: tidal velocity = v extracted from model at tide gauge point and time stamp in UTC
    tide_gauge = {}
    tide_gauge['fn'] = home+'HC_dolph_tidal_current_ATTideGauge_9445133.p'
    tide_gauge['dict'] = pickle.load(open(tide_gauge['fn'],'rb'))
    tide_gauge['df'] = pd.DataFrame({'dt_list':[dt for dt in tide_gauge['dict']['dt_list']]})
    tide_gauge['df']['datetime'] = pd.to_datetime(tide_gauge['df']['dt_list'])
    tide_gauge['df']['datetime'] = tide_gauge['df']['datetime'].dt.tz_localize(tz='UTC')

    # tide data: tidal velocity = v extracted from model at tide gauge point and time stamp in UTC
    tide_dolphinpen = {}
    tide_dolphinpen['fn'] = home+'HC_dolph_tidal_current_ATBangorDolphinPen.p'
    tide_dolphinpen['dict'] = pickle.load(open(tide_dolphinpen['fn'],'rb'))
    tide_dolphinpen['df'] = pd.DataFrame({'dt_list':[dt for dt in tide_dolphinpen['dict']['dt_list']]})
    tide_dolphinpen['df']['datetime'] = pd.to_datetime(tide_dolphinpen['df']['dt_list'])
    tide_dolphinpen['df']['datetime'] = tide_dolphinpen['df']['datetime'].dt.tz_localize(tz='UTC')

# tide data: tide height (m) at tide gauge
tide_height = {}
tide_height['fn'] = home+'NOAA_MLLW_Bangor_9445133.txt'
tide_height['df'] = pd.read_csv(tide_height['fn'],sep='\s+',engine='python',skiprows=13)
tide_height['df']['datetime'] =  pd.to_datetime(tide_height['df']['Date']+' '+tide_height['df']['Time'])
tide_height['df']['datetime'] = tide_height['df']['datetime'].dt.tz_localize(tz='UTC')

# prepare data series
log_process_mean=process['df']['process_term'][:]
log_DNA_mean = np.mean(np.log(ESP['df']['PB_quantity']))
process['df']['data0'] = np.exp(log_process_mean+log_DNA_mean)

ESP_trip_avg['df']['data0'] = ESP_trip_avg['df']['PB_quantity_mean']

dolphin['df']['data0'] = dolphin['df']['ndolphin_net']

tide_gauge['df']['data0'] = tide_gauge['dict']['v'][:]

tide_dolphinpen['df']['data0'] = tide_dolphinpen['dict']['v'][:]

tide_height['df']['data0'] = tide_height['df']['Pred']

# create dDNA/dt data
dDNAdt = {}
dDNAdt['df'] = pd.DataFrame()
myDNA = 'process'
if myDNA=='ESP':
    pt = ESP_trip_avg['df']['data0'].values
    dDNA = np.diff(pt)
    dt = np.diff(ESP_trip_avg['df']['datetime'].view('int64').values//10**9)
    dDNAdt['df']['data0'] = dDNA/dt
    pdt = ESP_trip_avg['df']['datetime'].values
    dDNAdt['df']['datetime'] = [pdt[i]+np.timedelta64(int(dt[i]),'s') for i in range(len(pdt)-1)]
else:
    pt = process['df']['data0'].values
    dDNA = np.diff(pt)
    dt = np.diff(process['df']['datetime'].view('int64').values//10**9)
    dDNAdt['df']['data0'] = dDNA/dt
    pdt = process['df']['datetime'].values
    dDNAdt['df']['datetime'] = [pdt[i]+np.timedelta64(int(dt[i]),'s') for i in range(len(pdt)-1)]


#~~~~~~~~~~~~~~~~~~~~

dataset_list = ESP_trip_avg,dolphin,tide_dolphinpen,tide_gauge,dDNAdt,tide_height,process

# generate master list bookends while converting to timestamps
start_datetime_list = []
end_datetime_list = []
for dataset in dataset_list:
    df = dataset['df']
    start_datetime_list.append(df.datetime.dt.floor('H').values[0])
    end_datetime_list.append(df.datetime.dt.ceil('H').values[-1])
    df['timestamp'] = df.datetime.values.view('int64') // 10 ** 9 # convert to time stamp for easier comparisons and interpolation

# generate master time list
master = pd.DataFrame()
start_datetime = max(start_datetime_list)
end_datetime = min(end_datetime_list)
master['Datetime'] = pd.date_range(start_datetime,end_datetime,freq='1H')
master['Timestamp (sec)'] = master['Datetime'].view('int64') // 10 ** 9

#~~~~~~~~~~~~~~~~~~~~~

# interpolate data onto master timestamps

plotting=False
if plotting:
    fig = plt.figure(figsize=(8,8))
    ax = fig.gca()
for dataset in dataset_list:
    # note: not saving in dataframes bc it will have different time indexes!
    dataset['data'] = np.interp(np.array(master['Timestamp (sec)']),np.array(dataset['df']['timestamp']),dataset['df']['data0'])
    
    if plotting:
        ax.plot(master['Datetime'],dataset['data'])

if plotting:
    # ax.scatter(ESP['df']['datetime'],ESP['df']['PB_quantity'])

    ax.set_yscale('log')
    plt.show(block=False)
    plt.pause(0.1)

else:
    master['DNA (copies/uL) (ESP triplicate mean)'] = ESP_trip_avg['data']
    master['N (dolphins)'] = dolphin['data']
    master['Current velocity at tide gauge (m/s)'] = tide_gauge['data']
    master['Current velocity at dolphin pen (m/s)'] = tide_dolphinpen['data']
    master['dDNAdt (copies/uL/s)'] = dDNAdt['data']
    master['Tide height (m)'] = tide_height['data']
    master['smoothed DNA (copies/uL) (observation error removed)'] = process['data']


    outfn = home+'ESP_box_model_terms_tidalheight_{}.csv'.format(tidecurrent)
    master.to_csv(outfn)