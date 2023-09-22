import os
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from PIL import Image
import numpy as np
from datetime import datetime, timedelta
import pytz
import netCDF4 as nc
import pandas as pd

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)

plt.close('all')
home = '/data2/pmr4/eab32/'
figname = home+'LO_data/eDNA/bangor_gmaps.png'

#find June sampling particle tracks
track_dir0 = home+'LO_output/tracks/'

f_list = os.listdir(track_dir0)
f_list.sort()
f_list = [x for x in f_list if (x[:11]=='hc_dolph_3d')&(x[-4]=='6')]

#set reference time
Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

first_sample_utc = Pacific.localize(datetime(2023,6,5,10,30)).astimezone(utc)
first_release = first_sample_utc-timedelta(days=2)
last_sample_utc = Pacific.localize(datetime(2023,6,5,15,30)).astimezone(utc)
dt_list0 = pd.date_range(start=first_release,end = last_sample_utc, freq="15min").to_pydatetime().tolist()

figcount =0
for dt in dt_list0:
    fig = plt.figure(figsize=(8, 8))

    # this is from the cartopy docs
    img_extent = (-122.754167, -122.695833, 47.718611, 47.788056)
    ax = plt.axes(projection=ccrs.PlateCarree())

    with Image.open(figname) as img:
        ax.imshow(img, extent=img_extent,transform=ccrs.PlateCarree())
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        gl=ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)
        gl.top_labels = False
        gl.right_labels = False

        for f in f_list:
    
            track_dir = track_dir0+f
    
            # build a keyname from the release filename
            file_list = os.listdir(track_dir)
            file_list = [x for x in file_list if x[:3]=='rel']
            rel_fn = file_list[0]
    
            ds = nc.Dataset(track_dir+'/'+rel_fn)
    
            ot = ds['ot'][:].data
    
            ot0 = utc.localize(datetime(1970,1,1)+timedelta(seconds=ot[0]))
            ot1 = utc.localize(datetime(1970,1,1)+timedelta(seconds=ot[-1]))
    
            if (dt>ot0)&(dt<ot1):
                dt_list = [utc.localize(datetime(1970,1,1)+timedelta(seconds=ott)) for ott in ot]
                t = argnearest(dt_list,dt)
                zmask = ds['z'][t,:]>(ds['zeta'][t,:]-2)
                if np.sum(zmask)>0:
                    ax.scatter(ds['lon'][t,zmask],ds['lat'][t,zmask],c='w',s=1,alpha=0.05)
            ds.close()
        ax.axis(img_extent)
        ax.text(0.9,0.1,dt.astimezone(Pacific).strftime("%m/%d/%Y, %H:%M:%S PDT"),transform=ax.transAxes,ha='right',va='center',color='white')
        outfn = home+'etools/plots/June5_samples/figure_{}.png'.format(str(figcount).zfill(3))
        plt.savefig(outfn)
        print(f'Saved to {outfn}')
        plt.close()
        figcount+=1


