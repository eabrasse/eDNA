import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from PIL import Image
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.dates as mdates
import pytz
from datetime import datetime, timedelta

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'
tab10 = plt.get_cmap('tab10',10)

#import background map
figname = home+'data/bangor_gmaps.png'

#import tidal height
mllw_data_fn = '../data/9445133_MLLW_20230603-20230605.txt'
df = pd.read_csv(mllw_data_fn,engine='python',skiprows=13,delimiter=r'\t+')
df['datetime']=pd.to_datetime(df['Date ']+' '+df.Time)

fig = plt.figure(figsize=(12, 8))
gs = GridSpec(3,2)

# this is from the cartopy docs
img_extent = (-122.754167, -122.695833, 47.718611, 47.788056)
ax = plt.subplot(gs[:,0],projection=ccrs.PlateCarree())

with Image.open(figname) as img:
    ax.imshow(img, extent=img_extent,transform=ccrs.PlateCarree())
    ax.plot([-122.74,-122.73],[47.73,47.74])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    gl=ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False

# add tide plot
axt = plt.subplot(gs[0,1])
axt.plot(df.datetime,df.Pred-np.mean(df.Pred),color=tab10(0))
# axt.set_xlabel('Time (UTC)')
axt.set_ylabel('Tide elevation (m)')
axt.set_xticklabels([''])
# axt.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d %H:%m"))
# plt.setp( axt.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
axt.grid()
ylim=axt.get_ylim()
xlim=axt.get_xlim()
axt.text(0.1,0.9,'a) Tidal height',transform=axt.transAxes)

#set reference time
Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc
first_sample_utc = Pacific.localize(datetime(2023,6,5,10,30)).astimezone(utc)
first_release = first_sample_utc-timedelta(days=2)
last_sample_utc = Pacific.localize(datetime(2023,6,5,15,30)).astimezone(utc)

axt.fill_between([first_sample_utc,last_sample_utc],[ylim[0]-1,ylim[0]-1],[ylim[1]+1,ylim[1]+1],color=tab10(0),alpha=0.25)
axt.set_ylim(ylim)

# add number of particles released
axp = plt.subplot(gs[1,1])
dt_list0 = pd.date_range(start=first_release,end = last_sample_utc, freq="15min").to_pydatetime().tolist()
particle_count = [100000*(1+(dt0-dt_list0[0]).days*24+divmod((dt0-dt_list0[0]).seconds,3600)[0]) for dt0 in dt_list0]
axp.plot(dt_list0,particle_count,color=tab10(8))
axp.set_xlim(xlim)
axp.set_ylabel('Particles')
axp.set_xlabel('Time (UTC)')
axp.text(0.1,0.9,'b) Number of particles released',transform=axp.transAxes)
axp.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d %H:%m"))
plt.setp( axp.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
axp.grid()

plt.show(block=False)
plt.pause(0.1)

