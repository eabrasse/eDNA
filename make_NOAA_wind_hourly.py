# -*- coding: utf-8 -*-
"""
modified from Parker's code
(c) Elizabeth Brasseale 11/21/2022

"""
import os
import sys
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
import csv

# set up filepaths
dir0='/Users/elizabethbrasseale/Projects/eDNA/'

obs_fn = dir0+'data/IOOS_Wind_Bremerton9445958_20241015_20241018.csv'
dfo = pd.read_csv(obs_fn,sep=',',engine='python',skiprows=[1])
dfo['dt'] = pd.to_datetime(dfo['time'])
dfo['Wind_Direction (rad)'] = dfo.apply(lambda row: ((270-row['Wind_Direction'])*np.pi/180),axis=1)
dfo['Uwind'] = dfo.apply(lambda row: (row['Wind_Speed']*np.cos(row['Wind_Direction (rad)'])),axis=1)
dfo['Vwind'] = dfo.apply(lambda row: (row['Wind_Speed']*np.sin(row['Wind_Direction (rad)'])),axis=1)

dto_list = pd.date_range(start=dfo['dt'].min().round('60min'),end=dfo['dt'].max().round('60min'),freq='1h',inclusive='both')
Uwind_1h = np.zeros(len(dto_list))
Vwind_1h = np.zeros(len(dto_list))
t=0
for dt in dto_list:
    dt0 = dt - timedelta(minutes=30)
    dt1 = dt + timedelta(minutes=30)
    
    gg = dfo[(dfo['dt']>=dt0)&(dfo['dt']<dt1)]
    if len(gg)<9: # 1 missing data point isn't a big deal. >1 missing point might be...
        Uwind_1h[t]=np.nan
        Vwind_1h[t]=np.nan
    else:
        Uwind_1h[t]=np.mean(gg['Uwind'])
        Vwind_1h[t]=np.mean(gg['Vwind'])
    t+=1

# save hourly mean wind to csv
hourlydict = {'dt_list':dto_list, 'Uwind':Uwind_1h, 'Vwind':Vwind_1h}
dfo_new = pd.DataFrame.from_dict(hourlydict)
dfo_new = dfo_new.dropna()
outfn = dir0+'data/IOOS_Wind_Bremerton9445958_20241015_20241018_hourly.csv'
dfo_new.to_csv(outfn,sep=',')
