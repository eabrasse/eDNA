#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot cutthroat trout DNA concentrations
"""

# importing modules, the beginning of all python code
import os
import sys
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import pandas as pd
import string
import numpy as np
import pytz

# some helpful plotting commands
plt.close('all')
tab10 = plt.get_cmap('tab10',10)
rainbow = plt.get_cmap('rainbow')

# load in data
data_fn = '/Users/elizabethbrasseale/Projects/eDNA/data/MURI_Module1_datasheets - 03_ESP1_Feb2023.csv'

# generate figure and axis handles
# fig,axs = plt.subplots(nrows=3, ncols=1,figsize=(8,6))
fig = plt.figure(figsize=(8,6))
ax = plt.gca()

df = pd.read_csv(data_fn,sep=',',engine='python')
df['date']=pd.to_datetime(df.ESP_date+' '+df.ESP_filter_t0)
df= df[df.date.notnull()]
df.date=df.date.dt.tz_localize(tz='America/Vancouver')
df.date=df.date.dt.tz_convert('UTC')

ax.scatter(df.date,df.PB_quantity,edgecolors=tab10(0),c='None')
udates = df.sort_values(by='date').date.unique()
PB_quantity_mean = np.zeros(udates.shape)
udates_rounded_to_hour = []
print('Times converted from PDT to UTC:')
count=0
for ud in udates:
    gg = df[df.date==ud]
    PB_quantity_mean[count] = gg.PB_quantity.mean()
    ax.scatter(ud,PB_quantity_mean[count],c=tab10(0))
    udates_rounded_to_hour.append(ud.round('60min'))
    print(ud.round('60min'))
    count+=1
ax.set_yscale('log')
ax.grid()
ax.set_xlabel('Time (UTC)')
ax.set_ylabel(r'Conc (copies/$\mathrm{\mu L}$)')
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d %H:%m"))
plt.setp( ax.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
Tues2pm = datetime(2023,1,31,14,0,0)
pdt = pytz.timezone('America/Vancouver')
Tues2pm = pdt.localize(Tues2pm)
Tues2pm = Tues2pm.astimezone(pytz.utc)
ax.axvline(x=Tues2pm,color=tab10(1))
ax.text(Tues2pm,1e4,'Tues 2pm',color=tab10(1))

# save hourly data
D = {'date':udates_rounded_to_hour,'PB_quantity_mean':PB_quantity_mean}
df_hourly_means = pd.DataFrame(data=D)
outfn = '/Users/elizabethbrasseale/Projects/eDNA/data/03_ESP1_Feb2023_hourly.csv'
df_hourly_means.to_csv(outfn)

# show plot
plt.tight_layout()
plt.show(block=False)
plt.pause(0.1)
