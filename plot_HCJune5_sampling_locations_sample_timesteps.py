import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from PIL import Image
import numpy as np
from datetime import datetime, timedelta
import pytz
import netCDF4 as nc
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.dates as mdates

def argnearest(items, pivot):
    td = [np.abs(item-pivot) for item in items]
    return np.argmin(td)

plt.close('all')
home = '/data2/pmr4/eab32/'

tab10 = plt.get_cmap('tab10',10)
Paired = plt.get_cmap('Paired',12)
# custom 2-color colormap
nocol = Paired(4)
yescol = Paired(2)
yesno_cmap = matplotlib.colors.ListedColormap([nocol,yescol])
yesno_col_list = [Paired(5),Paired(3)]

# load background map
figname = home+'LO_data/eDNA/bangor_gmaps.png'
img_extent = (-122.754167, -122.695833, 47.718611, 47.788056)

#find June sampling particle tracks
track_dir0 = home+'LO_output/tracks/'

f_list = os.listdir(track_dir0)
f_list.sort()
f_list = [x for x in f_list if (x[:28]=='hc_bangor_perimeter_point_3d')&(x[-4]=='6')]

#set reference time
Pacific = pytz.timezone("PST8PDT")
utc = pytz.utc

first_sample_utc = Pacific.localize(datetime(2023,6,5,10,30)).astimezone(utc)
first_release = first_sample_utc-timedelta(days=2)
last_sample_utc = Pacific.localize(datetime(2023,6,5,15,30)).astimezone(utc)
dt_list0 = pd.date_range(start=first_sample_utc,end = last_sample_utc, freq="15min").to_pydatetime().tolist()

#import sampling locations
sampling_data_fn = home+'LO_data/eDNA/June5_sampling_locations.csv'
dfs = pd.read_csv(sampling_data_fn,sep=',',engine='python')
dfs['approx_dt']=[Pacific.localize(pd.to_datetime(dfs.date_collected.values[i]+' '+dfs.approx_time_collected.values[i],utc=False)).astimezone(utc) for i in range(len(dfs.date_collected.values))]
dfs['mid_dt']=[Pacific.localize(pd.to_datetime(dfs.date_collected.values[i]+' '+dfs.mid_time_collected.values[i],utc=False)).astimezone(utc) for i in range(len(dfs.date_collected.values))]

#import tidal height
mllw_data_fn = home+'LO_data/eDNA/9445133_MLLW_20230603-20230605.txt'
dft = pd.read_csv(mllw_data_fn,engine='python',skiprows=13,delimiter=r'\t+')
dft['datetime']=[utc.localize(pd.to_datetime(dft['Date '].values[i]+' '+dft.Time.values[i])) for i in range(len(dft.Time.values))]

# import tidal current
current_fn = home+'LO_data/eDNA/ks0201_currents_20230603-20230605.txt'
dfks = pd.read_csv(current_fn,engine='python',skiprows=9,delimiter=r'\s+',names = ['Date','Time','Speed (cm/s)','Dir (true)'])
dfks['datetime']=[utc.localize(pd.to_datetime(dfks.Date.values[i]+' '+dfks.Time.values[i])) for i in range(len(dfks.Time.values))]
lon_ks0201 = -122 - 43.8/60
lat_ks0201 = 47 + 45.5/60

figcount =0

with Image.open(figname) as img:
    
    for dt in dt_list0:
        print(f'Timestep {(figcount+1)} of {len(dt_list0)}')
        # dt = Pacific.localize(datetime(2023,6,5,11,45)).astimezone(utc)
        dtts = pd.Timestamp(dt) #useful for comparisons with pandas dataframes
    
        fig = plt.figure(figsize=(12, 8))
        gs = GridSpec(3,2)

        ax = plt.subplot(gs[:,0],projection=ccrs.PlateCarree())

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
                    ax.scatter(ds['lon'][t,zmask],ds['lat'][t,zmask],c='w',s=1,alpha=0.05,zorder=100)
            ds.close()
        ax.axis(img_extent)
        
        ax.text(0.9,0.1,dt.astimezone(Pacific).strftime("%m/%d/%Y, %H:%M:%S PDT")+'\narrow = current, circle = 50 cm/s',transform=ax.transAxes,ha='right',va='center',color='white')

        # add tide plot
        axt = plt.subplot(gs[0,1])
        axt.plot(dft.datetime,dft.Pred-np.mean(dft.Pred),color=tab10(0))
        t = argnearest(dft.datetime,dtts)
        axt.plot(dft.datetime.values[t],dft.Pred[t]-np.mean(dft.Pred),marker='o',mec='k',mfc=tab10(0),markersize=10)
        # axt.set_xlabel('Time (UTC)')
        axt.set_ylabel('Tide elevation (m)')
        axt.set_xticklabels([''])
        # axt.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d %H:%m"))
        # plt.setp( axt.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
        axt.grid()
        ylim=axt.get_ylim()
        # xlim=axt.get_xlim()
        axt.text(0.1,0.9,'a) Tidal height',transform=axt.transAxes)

        # axt.fill_between([first_sample_utc,last_sample_utc],[ylim[0]-1,ylim[0]-1],[ylim[1]+1,ylim[1]+1],color=tab10(0),alpha=0.25)
        axt.set_ylim(ylim)
    
        # add tidal current
        t = argnearest(dfks.datetime,dtts)
        u = dfks['Speed (cm/s)'][t]*np.cos(dfks['Dir (true)'][t]*np.pi/180)
        v = dfks['Speed (cm/s)'][t]*np.sin(dfks['Dir (true)'][t]*np.pi/180)
        rad = 0.005
        current_legend = matplotlib.patches.Circle((lon_ks0201,lat_ks0201),rad,ec='white',fc='None',linewidth=1)
        ax.add_patch(current_legend)
        ax.quiver(lon_ks0201,lat_ks0201,u,v,pivot='tail',scale=50/rad,scale_units='xy',color='white')
    
    
        # add sample locations to map and sampling plot
        axs = plt.subplot(gs[1,1])
        axs.scatter(dfs.long[dfs.approx_dt>dt],dfs.lat[dfs.approx_dt>dt],marker='o',s=100,edgecolors='white',facecolors='None',zorder=200)
        axs.plot(dfs.mid_dt[dfs.approx_dt>dt],dfs.eDNA_present[dfs.approx_dt>dt],linestyle='None',marker='o',mec='gray',mfc='None',markersize=10)
        if dt>=first_sample_utc:
            ax.scatter(dfs[dfs.approx_dt<dt]['long'],dfs[dfs.approx_dt<dt]['lat'],marker='o',s=100,edgecolors='k',c=dfs.eDNA_present[dfs.approx_dt<dt],cmap=yesno_cmap,zorder=200)
            axs.scatter(dfs.mid_dt[dfs.approx_dt<dt],dfs.eDNA_present[dfs.approx_dt<dt],marker='o',s=100,edgecolors='k',c=dfs.eDNA_present[dfs.approx_dt<dt],cmap=yesno_cmap,zorder=200)
            if np.any(dfs.approx_dt==dtts):
                t=argnearest(dfs.approx_dt,dtts)
                ax.plot(dfs.long[t],dfs.lat[t],marker='o',markersize=12,mec='k',mfc=yesno_col_list[int(dfs.eDNA_present[t])],zorder=201)
                axs.plot(dfs.mid_dt[t],dfs.eDNA_present[t],linestyle='None',marker='o',mfc=yesno_col_list[int(dfs.eDNA_present[t])],mec='k',markersize=12)
        # axs.set_xlim(xlim)
        axs.set_xlabel('Time (UTC)')
        axs.set_ylabel('')
        axs.set_yticks([0,1])
        axs.set_yticklabels(['No','Yes'])
        ylim = axs.get_ylim()
        axs.set_ylim([ylim[0]-0.2,ylim[1]+0.5])
        axs.text(0.1,0.9,'b) eDNA present at sample?', transform=axs.transAxes)
        axs.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%-d %H:%m"))
        plt.setp( axs.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
        axs.grid()
    
        for ax in axt, axs:
            ax.axvline(x=dt,linestyle='dashed',linewidth=1,color='k')
            ax.set_xlim([first_sample_utc,last_sample_utc])


        outfn = home+'etools/plots/June5_perimeter_point_sampling_locs_timesteps/figure_{}.png'.format(str(figcount).zfill(3))
        plt.savefig(outfn)
        print(f'Saved to {outfn}')
        plt.close()
        figcount+=1


