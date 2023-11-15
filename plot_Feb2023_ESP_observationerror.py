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
data_fn = '/Users/elizabethbrasseale/Projects/eDNA/data/'

ESP_fn = data_fn + 'MURI_Module1_datasheets - 03_ESP1_Feb2023.csv'
proc_fn = data_fn + 'process_term_mean.csv'

# generate figure and axis handles
# fig,axs = plt.subplots(nrows=3, ncols=1,figsize=(8,6))
fig = plt.figure(figsize=(8,6))
ax = plt.gca()

df = pd.read_csv(ESP_fn,sep=',',engine='python')
df['date']=pd.to_datetime(df.ESP_date+' '+df.ESP_filter_t0)
df= df[df.date.notnull()]
df.date=df.date.dt.tz_localize(tz='America/Vancouver')
df.date=df.date.dt.tz_convert('UTC')

ax.scatter(df.date,df.PB_quantity,edgecolors=tab10(0),c='None')

udates = df.dropna(subset=['ESP_date','PB_quantity']).sort_values(by='date').date.unique()
PB_quantity_mean = np.nanmean(np.log(df.PB_quantity))

dfp = pd.read_csv(proc_fn,sep=',',engine='python')
demeaned_log_proc = dfp.process_term
log_proc = demeaned_log_proc + PB_quantity_mean
proc = np.exp(log_proc)
ax.plot(udates,proc,marker='o',markeredgecolor='None',color=tab10(0))

ax.set_yscale('log')
ax.grid()
ax.set_xlabel('Time (UTC)')
ax.set_ylabel(r'DNA conc (copies/$\mathrm{\mu L}$)')
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d %H:%m"))
plt.setp( ax.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')


# show plot
plt.tight_layout()
plt.show(block=False)
plt.pause(0.1)
