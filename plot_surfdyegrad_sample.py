import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
# from datetime import datetime, timedelta
# import pytz
from matplotlib.gridspec import GridSpec
# import matplotlib.dates as mdates
import pickle
from matplotlib import ticker
import string
import efun
import cmocean as cmo
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# custom 2-color colormap
landcol = matplotlib.colors.to_rgba('lightgray')
seacol = (1.0,1.0,1.0,0)
cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

atoz = string.ascii_lowercase

plt.close('all')
home = '/Users/elizabethbrasseale/Projects/eDNA/'

data_fn = home+'data/HC_surfdye_0605.p'

#options and data
cutoff = 1e-3 #what the cutoff is
perimeter = 150 # how long of a transect to test
plotting=True # whether to output snapshot plots
t_range = range(18,24) # time steps to plot
nudge = 150 # plotting limits
plot_dyemax = False # whether to output plot of dyemax

#floating pen location
lon0 = -122.733651
lat0 = 47.740037

tab10 = plt.get_cmap('tab10',10)

Dp = pickle.load(open(data_fn,'rb'))

# gather some fields, for convenience
lonr = Dp['lonr'][:]
latr = Dp['latr'][:]
maskr = Dp['maskr'][:]

dye0 = np.ma.masked_where(Dp['dye']>10,Dp['dye'])
nt,ny,nx = np.shape(dye0)

xgrid,ygrid = efun.ll2xy(lonr,latr,lon0,lat0)

x00 = np.argmin(np.abs(xgrid[0,:]))
y00 = np.argmin(np.abs(ygrid[:,0]))

xcenter=-80
ycenter=-160

xci = np.argmin(np.abs(xgrid[0,:]-xcenter))
yci = np.argmin(np.abs(ygrid[:,0]-ycenter))

spacing = 100

xpos_list = [xcenter-spacing, xcenter, xcenter+spacing]
ypos_list = [[ycenter-spacing],[ycenter],[ycenter+spacing]]
ni = len(xpos_list)
nj = len(ypos_list)

xpos_array = np.tile(np.array(xpos_list),(nj,1))
xinds_list = [np.argmin(np.abs(xgrid[0,:]-xpos)) for xpos in xpos_list]
xinds_array = np.tile(np.array(xinds_list),(nj,1)) #increase w/ column, same w/ row
xcount_array = np.tile(np.array(range(ni)),(nj,1))


ypos_array = np.tile(np.array([[ypos] for ypos in ypos_list]),(1,ni))
yinds_list = [[np.argmin(np.abs(ygrid[:,0]-ypos))] for ypos in ypos_list] #increase w/ row, same w/ column
yinds_array = np.tile(np.array(yinds_list),(1,ni))
ycount_array = np.tile(np.array([[j] for j in range(nj)]),(1,ni))


dye = np.zeros((ni,nj,nt))
mean_dye = np.zeros((ni,nj))

fw,fh = efun.gen_plot_props()

fig1 = plt.figure(figsize=(1.2*fw,fh))
# gs = GridSpec(nj,int(2*ni))
# fig2 = plt.figure(figsize=(fw*12,fh*0.5))

axmap = fig1.gca()

axmap.pcolormesh(xgrid,ygrid,maskr,cmap=cmap_mask,shading='nearest')
axmap.plot([0],[0],marker='o',linestyle='none',color='r',zorder=5000)
# ylim_list = []
axs = []

# for station in range(nstations):
for i in range(ni):
    ii = xinds_list[i]
    for j in range(nj):    
        axmap.plot(xpos_list[i],ypos_list[j],marker='o',linestyle='None',mec='k',mfc='white',zorder=50)

        jj =yinds_list[j]
    
        dye[i,j,:] = dye0[:,jj,ii].squeeze()
        mean_dye[i,j] = dye0[:,jj,ii].mean()
    

    
        # ax = fig1.add_subplot(gs[-(j+1),(i+ni)])
        # ax.plot(range(nt),dye[i,j,:],color='gray',lw=1.0)
        # ax.axhline(y=mean_dye[i,j],lw=1.5,linestyle='dashed',color='k')
        # ylim_list.append(ax.get_ylim())
        # ax.text(0.1,0.9,f'x={xpos_list[i]} m, y={ypos_list[j][0]}',transform=ax.transAxes,ha='left',va='top')
        # axs.append(ax)
        # ax.set_yscale('log')
        # if i>0:
        #     # ax.get_yaxis().set_visible(False)
        #     ax.set_yticklabels([])
        # if j>0:
        #     # ax.get_xaxis().set_visible(False)
        #     ax.set_xticklabels([])

ddyedx = np.zeros((ni,nj))
ddyedy = np.zeros((ni,nj))
for i in range(ni):
    for j in range(nj):
        if i==0:
            ddyedx[i,j] = (mean_dye[i+1,j]-mean_dye[i,j])/spacing
        elif i==(ni-1):
            ddyedx[i,j] = (mean_dye[i,j]-mean_dye[i-1,j])/spacing
        else:
            ddyedx[i,j] = (mean_dye[i+1,j]+mean_dye[i-1,j]-2*mean_dye[i,j])/(2*spacing)
    
        if j==0:
            ddyedy[i,j] = (mean_dye[i,j+1]-mean_dye[i,j])/spacing
        elif j==(nj-1):
            ddyedy[i,j] = (mean_dye[i,j]-mean_dye[i,j-1])/spacing
        else:
            ddyedy[i,j] = (mean_dye[i,j+1]+mean_dye[i,j-1]-2*mean_dye[i,j])/(2*spacing)
            
axmap.quiver(xpos_array.flatten(),ypos_array.flatten(),ddyedx.flatten(),ddyedy.flatten(),zorder=500)
axmap.set_xlabel('Dist from floating pen (m)')
axmap.set_ylabel('Dist from floating pen (m)')
# yl0 = np.min(ylim_list)
# yl1 = np.max(ylim_list)
# for ax in axs:
#     ax.set_ylim([1e-8,yl1])
#     ax.grid()

axmap.axis([-300,300,-300,300])
axmap.set_aspect(1)

fig1.subplots_adjust(right=0.9,left=0.2,bottom=0.05,top = 0.98,wspace=0.5,hspace=0.2)

plt.show(block=False)
plt.pause(0.1)
# plt.close()


