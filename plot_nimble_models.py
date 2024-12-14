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
# list all files in results folder
f_list = os.listdir(model_dir)
f_list.sort()
f_list = [x for x in f_list if (x[:2]=='m3')] # use only m3 run results

# initialize figure
fw,fh = efun.gen_plot_props()
fig = plt.figure(figsize=(2.5*fw,1.5*fh))
gs = GridSpec(3,3)

D = {}

# relate nimble model results filename, nimble model run title, axis to plot, and color in a dictionary
D['m30lnormobservation_plotting.csv'] = {'label':'M1 Constant Production\n and Loss','gs_pos':(0,0),'col':tab10(0)}
D['m31lineardolphin_plotting.csv'] = {'label':'M2 Linear dolphin','gs_pos':(1,0),'col':tab10(1)}
D['m32expdolphin_plotting.csv'] = {'label':'M3 Exponential\ndolphin','gs_pos':(1,1),'col':tab10(1)}
D['m33mixeddolphin_plotting.csv'] = {'label':'M4 Mixed-state\ndolphin','gs_pos':(1,2),'col':tab10(1)}
D['m34tidalloss_plotting.csv'] = {'label':'M5 Tidal loss','gs_pos':(2,0),'col':tab10(2)}
D['m35randomloss_plotting.csv'] = {'label':'M6 Random loss','gs_pos':(2,1),'col':tab10(2)}

# open the ESP dataset used to fit the model
ESPfn = home+'code/R/nimble/Data/ESP_timestamps_mLseawater.csv'
ESP = pd.read_csv(ESPfn,sep=',')
ts0 = ESP['timestamp_0'].min()
xscale = 1/3600

count=0
# loop through nimble models and plot each against ESP data on the axis defined above
for fname in f_list:
    #read in nimble model results
    df = pd.read_csv(model_dir+fname,sep=',',engine='python')
    fdict = D[fname]
    ax = fig.add_subplot(gs[fdict['gs_pos']]) #grab the axis
    
    # because of log axis, standard deviation is asymmetric. here we just plot upper limit
    yplus = df.model+df.sigma

    
    # plot ESP data using bars with width = sampling integration time as in previous figure
    for i in ESP.index:
        ax.fill_between([xscale*(ESP['timestamp_0'][i]-ts0),xscale*(ESP['timestamp_1'][i]-ts0)],[ESP['copies_per_mLseawater'][i]]*2,
            alpha=1.0,color='k',zorder=2,edgecolor='none')#,lw=4)
    
    # then plot the nimble model output as a dashed line plot over the bars
    # note that conceptual model should represent continuous time series, and ESP should subselect
    ax.plot(df.hours_since_ESP_start,df.model,color=fdict['col'],marker='None',linestyle='dashed',linewidth=2,zorder=150)
    # fill to indicate uncertainty
    ax.fill_between(df.hours_since_ESP_start,yplus,color=fdict['col'],alpha=0.15,ec='none')
    
    # set scale to log
    ax.set_yscale('log')
    
    # add label
    ax.text(0.1,0.9,fdict['label'],transform=ax.transAxes,color=fdict['col'],fontsize=10,va='top',ha='left',fontweight='bold',zorder=300)
    
    # add y-label to only leftmost plots to keep things uncluttered
    if fdict['gs_pos'][1]==0:
        ax.set_ylabel('DNA conc\n'+r'(copies $\mathrm{mL}^{-1}$)')
    else:
        ax.set_yticklabels(['' for yt in ax.get_yticks()])
    
    # add x-label to only bottom-most plots to keep things uncluttered
    if (fdict['gs_pos'][0]==2) or (fdict['gs_pos'][1]==2):
        ax.set_xlabel('Time since sampling start (hr)')
    else:
        ax.set_xticklabels(['' for xt in ax.get_xticks()])
    
    #grid helps clarity
    ax.grid(color='lightgray')
    
    # add a legend once
    if count==0:
        ax.text(0.1,0.7,'Observations',transform=ax.transAxes,color='k',fontsize=9,fontweight='normal',ha='left',va='top')
    
    count+=1
    

fig.subplots_adjust(top=0.98,left=0.1,right=0.98,bottom=0.1,hspace=0.1,wspace=0.1)
# plt.show(block=False)
# plt.pause(0.1)
outfn = home+'plots/nimble_models_mLseawater.png'
plt.savefig(outfn,format='png',dpi=400)
