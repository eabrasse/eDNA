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


###
# Pedro's mifish haplotype data

mifish = {}
mifish['fn'] = home+'ESP_ASVs_filtered.csv'
mifish['df'] = pd.read_csv(mifish['fn'],sep=',',engine='python')
mifish['df']['Sample code2'] = mifish['df']['Sample code'].str.slice(0,-5)
sample_list = [f'MFU.52345.{count:}' for count in range(171,213)]
mifish['df']=mifish['df'][mifish['df']['Sample code2'].isin(sample_list)]
column_names = list(mifish['df'].columns)
column_sub = [column for column in column_names if column not in ['Sample code','Sample code2','timestamp_0','timestamp_1']]
ndolphhap = len(column_sub)
for column in column_sub:
    mifish['df'][column] = mifish['df'][column].astype(float)
###
# count total reads in two ways

# first sum over mifish by row
mifish['df']['Nreads'] = mifish['df'].apply(lambda row: np.sum(np.array([row[column] for column in column_sub])),axis=1)

#then look at Eily's .csv
metamifish = {}
metamifish['fn'] = home+'ASV_table.csv'
metamifish['df'] = pd.read_csv(metamifish['fn'],sep=',',engine='python')
metamifish['df'] = metamifish['df'].groupby(metamifish['df'].Sample).sum()
sample_list = [f'MFU-52345-{count:}' for count in range(171,213)]
metamifish['df'] = metamifish['df'][metamifish['df'].index.isin(sample_list)]

###

fw,fh=efun.gen_plot_props()
fig = plt.figure(figsize=(fw*2,fh))
ax = fig.gca()

ax.scatter(range(171,213),mifish['df']['Nreads'],c=tab10(0),label='Pedro datasheet')
ax.scatter(range(171,213),metamifish['df']['nReads'],c=tab10(1),label='Eily datasheet')
ax.legend()
ax.set_xlabel('Sample')
ax.set_ylabel('N reads')

fig.subplots_adjust(bottom=0.2)
fig.show()
plt.pause(0.1)

