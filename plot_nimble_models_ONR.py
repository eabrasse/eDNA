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
import matplotlib.dates as mdates
from matplotlib.dates import DayLocator, HourLocator, MonthLocator,DateFormatter, drange
import string
import os

plt.close('all')
tab10 = plt.get_cmap('tab10',10)
tab20b = plt.get_cmap('tab20b',20)
# tz = pytz.utc
# pdt = pytz.timezone('America/Vancouver')
atoz = string.ascii_lowercase

# load in data
home = '/Users/elizabethbrasseale/Projects/eDNA/'
model_dir = home+'code/R/nimble/Results/'

f_list = os.listdir(model_dir)
f_list.sort()
# f_list = [x for x in f_list if (x[:2]=='m3')]
f_list = ['m33mixeddolphin_plotting.csv']


fw,fh = efun.gen_plot_props()
fig = plt.figure(figsize=(1.5*fw,0.8*fh))
# gs = GridSpec(3,3)

D = {}

# add properties
# D['m30lnormobservation_plotting.csv'] = {'label':'M1 Constant Production\n and Loss','gs_pos':(0,0),'col':tab10(0)}
# D['m31lineardolphin_plotting.csv'] = {'label':'M2 Linear dolphin','gs_pos':(1,0),'col':tab10(1)}
# D['m32expdolphin_plotting.csv'] = {'label':'M3 Exponential\ndolphin','gs_pos':(1,1),'col':tab10(1)}
D['m33mixeddolphin_plotting.csv'] = {'label':'Mixed-state\neDNA production model','gs_pos':(1,2),'col':tab10(1)}
# D['m34tidalloss_plotting.csv'] = {'label':'M5 Tidal loss','gs_pos':(2,0),'col':tab10(2)}
# D['m35randomloss_plotting.csv'] = {'label':'M6 Random loss','gs_pos':(2,1),'col':tab10(2)}

ESPfn = home+'code/R/nimble/Data/ESP_timestamps_mLseawater.csv'
ESP = pd.read_csv(ESPfn,sep=',')
ts0 = ESP['timestamp_0'].min()
xscale = 1/3600

count=0
for fname in f_list:
    df = pd.read_csv(model_dir+fname,sep=',',engine='python')
    fdict = D[fname]
    # ax = fig.add_subplot(gs[fdict['gs_pos']])
    ax = fig.gca()
    
    yplus = df.model+df.sigma
    # ymin = df.model-df.sigma
    # ymin[ymin<0] = 1
    
    # ax.plot(df.hours_since_ESP_start,df.ESP_DNA_conc,marker='None',mfc='none',mec='none',markersize=5,linestyle='solid',linewidth=1,color='k',zorder=100)
    for i in ESP.index:
        ax.fill_between([xscale*(ESP['timestamp_0'][i]-ts0),xscale*(ESP['timestamp_1'][i]-ts0)],[ESP['copies_per_mLseawater'][i]]*2,
            alpha=1.0,color='k',zorder=2,edgecolor='none')#,lw=4)

    ax.plot(df.hours_since_ESP_start,df.model,color=fdict['col'],marker='None',linestyle='dashed',linewidth=2,zorder=150)
    ax.fill_between(df.hours_since_ESP_start,yplus,color=fdict['col'],alpha=0.15,ec='none')
    
    ax.set_yscale('log')
    
    ax.text(0.1,0.9,fdict['label'],transform=ax.transAxes,color=fdict['col'],fontsize=10,va='top',ha='left',fontweight='bold',zorder=300)
    
    # if fdict['gs_pos'][1]==0:
    ax.set_ylabel('DNA conc\n'+r'(copies $\mathrm{mL}^{-1}$)')
    # else:
        # ax.set_yticklabels(['' for yt in ax.get_yticks()])
    
    
    # if (fdict['gs_pos'][0]==2) or (fdict['gs_pos'][1]==2):
    ax.set_xlabel('Time since sampling start (hr)')
    # else:
        # ax.set_xticklabels(['' for xt in ax.get_xticks()])
    ax.grid(color='lightgray')
    
    # if count==0:
    ax.text(0.1,0.7,'Observations',transform=ax.transAxes,color='k',fontsize=9,fontweight='normal',ha='left',va='top')
    
    count+=1
    

fig.subplots_adjust(top=0.98,left=0.2,right=0.98,bottom=0.15)
# plt.show(block=False)
# plt.pause(0.1)
outfn = home+'plots/nimble_models_mLseawater_ONR.png'
plt.savefig(outfn,format='png',dpi=400)
