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
# from scipy.stats import linregress

plt.close('all')
tab10 = plt.get_cmap('tab10',10)
tz = pytz.utc
pdt = pytz.timezone('America/Vancouver')

# load in data
home = '/Users/elizabethbrasseale/Projects/eDNA/data/'

datafn = home+'ESP_box_model_terms_tidalheightfrac_overheight.csv'
df = pd.read_csv(datafn,sep=',',engine='python')
df['Datetime']=pd.to_datetime(df['Datetime'])

height = 0.5*((1/df['Over height (1/m)'][:-1])+(1/df['Over height (1/m)'][1:]))

# lhs = height*df['DNA (copies/uL) (ESP triplicate mean)'][1:]-df['DNA (copies/uL) (ESP triplicate mean)'][:-1]*(0.23-0.5*(df['Tide height fractional change (unitless)'][1:]+df['Tide height fractional change (unitless)'][:-1]))
lhs = df['DNA (copies/uL) (ESP triplicate mean)'][:]
# ndolphins = 0.5*(df['N (dolphins)'][:-1]+df['N (dolphins)'][1:])
ndolphins = df['N (dolphins)'][:]

fw,fh = efun.gen_plot_props()
fig = plt.figure(figsize=(6,4))
ax = fig.gca()

ax.plot(ndolphins,lhs,linestyle='None',marker='o',mec='k',mfc = 'gray',markersize=5)

ax.set_xlabel(r'N dolphins')
# ax.set_yscale('log')
# ax.set_ylabel(r' $(H+SSH_{t-1})\times(\mu_{t}-\mu_{t-1}(\alpha - \frac{SSH_{t}-SSH_{t-1}}{H+SSH_{t-1}}))$')
ax.set_ylabel('DNA copies/uL')

plt.subplots_adjust(left=0.2,bottom=0.23,right=0.99,top=0.99)
plt.show(block=False)
plt.pause(0.1)
