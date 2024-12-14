import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
# import pytz
from matplotlib.gridspec import GridSpec
# import matplotlib.dates as mdates
import netCDF4 as nc
import pickle
from matplotlib import ticker
import string
import efun
import cmocean as cmo
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

tab10 = plt.get_cmap('tab10',10)

atoz = string.ascii_lowercase

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'

data_fn = home+'data/f2023.06.05_test_case/'
fn_list = os.listdir(data_fn)
fn_list.sort()
fnh_list = [x for x in fn_list if x[6:9]=='his']
hrs = len(fnh_list)

#floating pen location
floatingpen = {}
floatingpen['lon0'] = -122.733651
floatingpen['lat0'] = 47.740037

tidegauge = {}
tidegauge['lon0']=-(122 + 43.8/60)
tidegauge['lat0']=47 + 45.5/60

location_list = floatingpen, tidegauge
nlocs = len(location_list)

dt_list = []
tt=0
# for hr in range(hrs):
hr=0

oceanfnh = data_fn + fnh_list[hr]

dsh = nc.Dataset(oceanfnh)

lonr = dsh['lon_rho'][:]
latr = dsh['lat_rho'][:]
maskr = dsh['mask_rho'][:]
lonv = dsh['lon_v'][:]
latv = dsh['lat_v'][:]
maskv = dsh['mask_v'][:]
# lonvecr = lonr[0,:]
# latvecr = latr[:,0]
# i_r = efun.find_nearest_ind(lonvecr,floatingpen['lon0'])
# j_r = efun.find_nearest_ind(latvecr,floatingpen['lat0'])
v = dsh['v'][:]
vbar = dsh['vbar'][:]

h = dsh['h'][:]
fudge = 0.05

fig = plt.figure(figsize=(10,10))
ax = fig.gca()

p=ax.pcolormesh(lonr,latr,h,cmap='cividis_r',vmin=0,vmax=20,alpha=0.5)
ax.pcolormesh(lonv,latv,vbar[0,:],cmap=cmo.cm.balance,vmin=-0.5,vmax=0.5)
ax.axis([floatingpen['lon0']-fudge, floatingpen['lon0']+fudge, floatingpen['lat0']-fudge, floatingpen['lat0']+fudge])
contourlist = np.arange(0,50,2)
ncon = len(contourlist)
colors = plt.cm.cividis_r(np.linspace(0,1,ncon))
cs=ax.contour(lonr,latr,h,levels=contourlist,colors=colors,linewidths=[2])
ax.clabel(cs,inline=True,fontsize=10)
ax.plot(floatingpen['lon0'],floatingpen['lat0'],marker='*',mec='k',mfc='yellow',markersize=10)
cbar=plt.colorbar(p,ticks=contourlist)
cbar.add_lines(cs)
plt.show(block=False)
plt.pause(0.1)
#
# if tt==0:
#
#     #load in lat & lon on u & v grids
#     lonr = dsh['lon_rho'][:]
#     latr = dsh['lat_rho'][:]
#     maskr = dsh['mask_rho'][:]
#     lonvecr = lonr[0,:]
#     latvecr = latr[:,0]
#
#     lonu = dsh['lon_u'][:]
#     latu = dsh['lat_u'][:]
#     masku = dsh['mask_u'][:]
#     lonvecu = lonu[0,:]
#     latvecu = latu[:,0]
#
#     lonv = dsh['lon_v'][:]
#     latv = dsh['lat_v'][:]
#     maskv = dsh['mask_v'][:]
#     lonvecv = lonv[0,:]
#     latvecv = latv[:,0]
#
#     # S = zrfun.get_basic_info(oceanfnh,only_S=True)
#
#     for location in location_list:
# get tidal location index

#
#         location['i_u'] = efun.find_nearest_ind(lonvecu,location['lon0'])
#         location['j_u'] = efun.find_nearest_ind(latvecu,location['lat0'])
#
#         location['i_v'] = efun.find_nearest_ind(lonvecv,location['lon0'])
#         location['j_v'] = efun.find_nearest_ind(latvecv,location['lat0'])
#
#
#         location['h'] = dsh['h'][location['j_r'],location['i_r']]
#
#     # depth_list = [depth for depth in depth_list0 if np.abs(depth)<h]
#     # ndepth = len(depth_list)
#
#     #initialize arrays
#     zeta = np.zeros((hrs,nlocs))
#     usurf = np.zeros((hrs,nlocs))
#     ubar = np.zeros((hrs,nlocs))
#     vsurf = np.zeros((hrs,nlocs))
#     vbar = np.zeros((hrs,nlocs))
#
#     # if diag:
#     #     sustr = np.zeros((nt))
#     #     svstr = np.zeros((nt))
#
#
# ot = dsh['ocean_time'][:].data[0]
# dt_list.append(datetime(1970,1,1)+timedelta(seconds=ot))
#
# outfn0 = home+'plots/floating pen tide/tide_vbar{:}.png'.format(hr)
# fig=plt.figure(figsize=(8,8))
# ax = fig.gca()
# p=ax.pcolormesh(lonv,latv,dsh['vbar'][0,:,:],cmap=cmo.cm.balance,vmin=-0.5,vmax=0.5)
# ax.plot(floatingpen['lon0'],floatingpen['lat0'],'x',color=tab10(0))
# ax.plot(tidegauge['lon0'],tidegauge['lat0'],'x',color=tab10(1))
# ax.text(0.1,0.9,'VERTICALLY-AVERAGED N-S VELOCITY (m/s)',transform=ax.transAxes)
# plt.colorbar(p)
# plt.savefig(outfn0)
# plt.close()
#
# outfn0 = home+'plots/floating pen tide/tide_ubar{:}.png'.format(hr)
# fig=plt.figure(figsize=(8,8))
# ax = fig.gca()
# p=ax.pcolormesh(lonu,latu,dsh['ubar'][0,:,:],cmap=cmo.cm.balance,vmin=-0.5,vmax=0.5)
# ax.plot(floatingpen['lon0'],floatingpen['lat0'],'x',color=tab10(0))
# ax.plot(tidegauge['lon0'],tidegauge['lat0'],'x',color=tab10(1))
# ax.text(0.1,0.9,'VERTICALLY-AVERAGED E-W VELOCITY (m/s)',transform=ax.transAxes)
# plt.colorbar(p)
# plt.savefig(outfn0)
# plt.close()
#
# outfn0 = home+'plots/floating pen tide/tide_vsurf{:}.png'.format(hr)
# fig=plt.figure(figsize=(8,8))
# ax = fig.gca()
# p=ax.pcolormesh(lonv,latv,dsh['v'][0,-1,:,:],cmap=cmo.cm.balance,vmin=-0.5,vmax=0.5)
# ax.plot(floatingpen['lon0'],floatingpen['lat0'],'x',color=tab10(0))
# ax.plot(tidegauge['lon0'],tidegauge['lat0'],'x',color=tab10(1))
# ax.text(0.1,0.9,'SURFACE N-S VELOCITY (m/s)',transform=ax.transAxes)
# plt.colorbar(p)
# plt.savefig(outfn0)
# plt.close()
#
# outfn0 = home+'plots/floating pen tide/tide_usurf{:}.png'.format(hr)
# fig=plt.figure(figsize=(8,8))
# ax = fig.gca()
# p=ax.pcolormesh(lonu,latu,dsh['u'][0,-1,:,:],cmap=cmo.cm.balance,vmin=-0.5,vmax=0.5)
# ax.plot(floatingpen['lon0'],floatingpen['lat0'],'x',color=tab10(0))
# ax.plot(tidegauge['lon0'],tidegauge['lat0'],'x',color=tab10(1))
# ax.text(0.1,0.9,'SURFACE E-W VELOCITY (m/s)',transform=ax.transAxes)
# plt.colorbar(p)
# plt.savefig(outfn0)
# plt.close()
#
# for l in range(nlocs):
#     location=location_list[l]
#     zeta[tt,l] = dsh['zeta'][0,location['j_r'],location['i_r']]
#     ubar[tt,l] = dsh['ubar'][0,location['j_u'],location['i_u']]
#     usurf[tt,l] = dsh['u'][0,-1,location['j_u'],location['i_u']]
#     vbar[tt,l] = dsh['vbar'][0,location['j_v'],location['i_v']]
#     vsurf[tt,l] = dsh['v'][0,-1,location['j_v'],location['i_v']]
    # if diag:
    #     sustr[tt] = dsh['sustr'][0,j_u,i_u]
    #     svstr[tt] = dsh['svstr'][0,j_v,i_v]

# zvec = zrfun.get_z(np.array(h),np.array(zeta[tt]),S,only_rho=True)

# zcount=0
# for depth in depth_list:
#     k_r = zfun.find_nearest_ind(zvec,depth)
#     u[tt,zcount] = dsh['u'][0,k_r,j_u,i_u]
#     v[tt,zcount] = dsh['v'][0,k_r,j_v,i_v]
#     zcount+=1


dsh.close()
tt = tt+1
#
# # print('Finished {:} (day {:} of {:})'.format(dt_list[tt-1].strftime('%b %-d'),(t+1),ndays))
# label='zeta','ubar','u_surf','vbar','v_surf'
# count =0
# for var in [zeta,ubar,usurf,vbar,vsurf]:
#     fig=plt.figure(figsize=(8,8))
#     ax = fig.gca()
#     for l in range(nlocs):
#         ax.plot(dt_list,var[:,l],color=tab10(l))
#
#     ax.text(0.1,0.9,label[count],transform=ax.transAxes)
#     ax.text(0.1,0.8,'floating pen',transform=ax.transAxes,color=tab10(0))
#     ax.text(0.1,0.7,'tide gauge',transform=ax.transAxes,color=tab10(1))
#     ax.set_xlabel('date time')
#     outfn = home+'plots/floating pen tide/'+label[count]+'.png'
#
#     plt.savefig(outfn)
#     plt.close()
#     count+=1
# fig1.subplots_adjust(right=0.9,left=0.09,bottom=0.05,top = 0.98,wspace=0.5,hspace=0.2)

# plt.show(block=False)
# plt.pause(0.1)
# plt.close()


