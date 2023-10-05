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
ESP['df']['datetime'] = pd.to_datetime(ESP['df']['ESP_date']+' '+ESP['df']['ESP_filter_tf'])
ESP['df'] = ESP['df'].dropna(subset=['ESP_date','PB_quantity']) # ignore data without timestamps
process['df']['datetime']= ESP['df']['datetime'].unique() # remove replicates of timestamps
process['df']=process['df'].sort_values(by='datetime')
process['df']['datetime'] = process['df']['datetime'].dt.tz_localize(tz='America/Vancouver')
process['df']['datetime'] = process['df']['datetime'].dt.tz_convert('UTC')

# dolphin occupancy data, digitized from Navy clipboard
dolphin = {}
dolphin['fn'] = home+'dolphin_presence.csv'
dolphin['df'] = pd.read_csv(dolphin['fn'],sep=',',engine='python')
dolphin['df'] = dolphin['df'].dropna()
dolphin['df']['datetime'] = pd.to_datetime(dolphin['df'].date+' '+dolphin['df'].time)
dolphin['df']['datetime'] = dolphin['df']['datetime'].dt.tz_localize(tz='America/Vancouver')
dolphin['df']['datetime'] = dolphin['df']['datetime'].dt.tz_convert('UTC')

# tide data: tidal speed = sqrt(u^2+v^2) extracted from model at tide gauge point and time stamp in UTC
tide = {}
tide['fn'] = home+'HC_dolph_tidal_current_ATTideGauge_9445133.csv'
tide['df'] = pd.read_csv(tide['fn'],sep=',',engine='python')
tide['df']['datetime'] = pd.to_datetime(tide['df']['Date Time'])
tide['df']['datetime'] = tide['df']['datetime'].dt.tz_localize(tz='UTC')

# prepare data series
log_process_mean=process['df']['process_term'][:]
log_DNA_mean = np.mean(np.log(ESP['df']['PB_quantity']))
process['df']['data0'] = np.exp(log_process_mean+log_DNA_mean)

dolphin['df']['data0'] = dolphin['df']['ndolphin_net']

tide['df']['data0'] = tide['df']['vel (m/s)']

# create dDNA/dt data
dDNAdt = {}
dDNAdt['df'] = pd.DataFrame()
pt = process['df']['data0'].values
dDNA = np.diff(pt)
dt = np.diff(process['df']['datetime'].view('int64').values//10**9)
dDNAdt['df']['data0'] = dDNA/dt
pdt = process['df']['datetime'].values
dDNAdt['df']['datetime'] = [pdt[i]+np.timedelta64(int(dt[i]),'s') for i in range(len(pdt)-1)]


#~~~~~~~~~~~~~~~~~~~~

# generate master list bookends while converting to timestamps
start_datetime_list = []
end_datetime_list = []
for df in process['df'],dolphin['df'],tide['df'],dDNAdt['df']:
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

# fig = plt.figure(figsize=(8,8))
# ax = fig.gca()
for dataset in process,dolphin,tide,dDNAdt:
    # note: not saving in dataframes bc it will have different time indexes!
    dataset['data'] = np.interp(np.array(master['Timestamp (sec)']),np.array(dataset['df']['timestamp']),dataset['df']['data0'])
    
    # ax.plot(master_datetime,dataset['data'])

# ax.set_yscale('log')
# plt.show(block=False)
# plt.pause(0.1)

master['DNA (copies/uL)'] = process['data']
master['N (dolphins)'] = dolphin['data']
master['Tidal Current Speed (m/sec)'] = tide['data']
master['dDNAdt (copies/uL/sec)'] = dDNAdt['data']

outfn = home+'ESP_box_model_terms.csv'
master.to_csv(outfn)