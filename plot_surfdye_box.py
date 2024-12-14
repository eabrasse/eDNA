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

xgrid,ygrid = efun.ll2xy(lonr,latr,lon0,lat0)

x00 = np.argmin(np.abs(xgrid[0,:]))
y00 = np.argmin(np.abs(ygrid[:,0]))

upright = {'pax1':{'md_list':[],'color':tab10(0)},'pax2':{'md_list':[],'color':tab10(2)},'mfc':'m','lw':1,'ls':'dashed'}
diag = {'pax1':{'md_list':[],'color':tab10(0)},'pax2':{'md_list':[],'color':tab10(2)},'mfc':'b','lw':2,'ls':'solid'}


xp0 = np.argmin(np.abs(xgrid[0,:]+perimeter))
xp1 = np.argmin(np.abs(xgrid[0,:]-perimeter))
yp0 = np.argmin(np.abs(ygrid[:,0]+perimeter))
yp1 = np.argmin(np.abs(ygrid[:,0]-perimeter))

upright['pax1']['xinds'] = range(xp0,xp1+1)
upright['pax1']['yinds'] = [y00]*len(upright['pax1']['xinds'])

upright['pax2']['yinds'] = range(yp0,yp1+1)
upright['pax2']['xinds'] = [x00]*len(upright['pax2']['yinds'])

#this is how it should be...
xperimeter = perimeter*np.sqrt(2)*0.5
xpx0 = np.argmin(np.abs(xgrid[0,:]+xperimeter))
xpx1 = np.argmin(np.abs(xgrid[0,:]-xperimeter))
ypx0 = np.argmin(np.abs(ygrid[:,0]+xperimeter))
ypx1 = np.argmin(np.abs(ygrid[:,0]-xperimeter))

#conveniently, both xp and yp are 20 indices long
#exploit this for now
x_range = range(xpx0,xpx1+1)
y_range = range(ypx0,ypx1+1)

diag['pax1']['xinds'] = range(xpx0,xpx1+1)
diag['pax1']['yinds'] = range(ypx0,ypx1+1)

diag['pax2']['xinds'] = range(xpx0,xpx1+1)
diag['pax2']['yinds'] = range(ypx0,ypx1+1)[::-1]

for orientation in upright, diag:
    for pax in 'pax1','pax2':
        orientation[pax]['dist'] = np.array([np.sign(xgrid[orientation[pax]['yinds'][i],orientation[pax]['xinds'][i]])*np.sqrt(xgrid[orientation[pax]['yinds'][i],orientation[pax]['xinds'][i]]**2+ygrid[orientation[pax]['yinds'][i],orientation[pax]['xinds'][i]]**2) for i in range(len(orientation[pax]['xinds']))])
# distvec = np.array([np.sqrt(xgrid[0,x_range[i]]**2+ygrid[y_range[i],0]**2)*np.sign(xgrid[0,x_range[i]]) for i in range(len(x_range))])







fw,fh = efun.gen_plot_props()

if plotting:
    fig = plt.figure(figsize=(fw*3,fh*1.5))
    gs = GridSpec(3,len(t_range))




count=0
for t in t_range:
    
    dye = dye0[t,:]
    for orientation in upright, diag:
        for pax in 'pax1','pax2':
            orientation[pax]['dye'] = np.array([dye[orientation[pax]['yinds'][i],orientation[pax]['xinds'][i]] for i in range(len(orientation[pax]['xinds']))])

    
    #plot colormap of dye
    if plotting:
        axs = [fig.add_subplot(gs[0,count]), fig.add_subplot(gs[1,count]), fig.add_subplot(gs[2,count])]
        
        axs[0].pcolormesh(xgrid,ygrid,maskr,cmap=cmap_mask,shading='nearest',zorder=100)
        for orientation in upright, diag:
            orientation['pax1']['ax'] = axs[1]
            orientation['pax2']['ax'] = axs[2]
        p=axs[0].pcolormesh(xgrid,ygrid,dye,cmap=cmo.cm.matter,norm=matplotlib.colors.LogNorm(vmax=cutoff*1e2,vmin=cutoff*1e-2),zorder=1,shading='nearest')
    
        axs[0].contour(xgrid,ygrid,maskr,levels=[0.5],colors=['k'],linewidths=[1.5],zorder=150)
        pl=axs[0].contour(xgrid,ygrid,dye,colors=['k'],linewidths=[0.5],levels=[cutoff],zorder=2)
        # axs[0].set_xlabel('Dist from floating pen (m)')
        # axs[0].set_ylabel('Dist from floating pen (m)')
        if t == t_range[-1]:
            cbaxes = inset_axes(axs[0], width="4%", height="40%", loc='lower right',bbox_transform=axs[0].transAxes,bbox_to_anchor=(-0.28,0,1,1))
            cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
            cb.add_lines(pl)

    
        axs[0].set_aspect(1)
        axs[0].axis([-nudge,nudge,-nudge,nudge])
        axs[0].plot(0,0,marker='x',color='yellow')
        axs[0].text(0.1,0.9,'t = {} hr'.format(t),transform=axs[0].transAxes,ha='left',va='top',zorder=5000)

        for orientation in diag,upright:
            for pax in 'pax1','pax2':
                axs[0].plot(xgrid[orientation[pax]['yinds'],orientation[pax]['xinds']],ygrid[orientation[pax]['yinds'],orientation[pax]['xinds']],color=orientation[pax]['color'],lw=orientation['lw'],linestyle=orientation['ls'])
                orientation[pax]['ax'].plot(orientation[pax]['dist'],orientation[pax]['dye'],color=orientation[pax]['color'],lw=orientation['lw'],linestyle=orientation['ls'])
    
    
        axs[1].set_xticklabels([''])


        for ax in axs[1:]:
        
            ax.set_xticks([-100,-80,-60,-40,-20,20,40,60,80,100],minor=True)
            ax.grid(which='both')
            ax.set_yscale('log')
            # ax.set_ylim([0,0.03])
            ax.set_ylim([cutoff*1e-1,cutoff*1e2])
            ax.axhline(y=cutoff,linestyle='solid',color='k')
            if count>0:
                # ax.get_yaxis().set_visible(False)
                ax.set_yticklabels([''])
    
    
    for pa in 'pax1','pax2':
        for orientation in upright,diag:
            dyemask = [orientation[pa]['dye']>cutoff]
            md_ind = np.argmax(np.abs(orientation[pa]['dist'][dyemask])) # index of maximum distance away you can see dye > 1e-3
            md = orientation[pa]['dist'][dyemask][md_ind] # max distance, but with sign
            if plotting:
                orientation[pa]['ax'].plot(md,orientation[pa]['dye'][dyemask][md_ind],marker='*',mec='k',mfc=orientation['mfc'],markersize=10)
            orientation[pa]['md_list'].append(np.abs(md))

    
    count+=1
    

if plot_dyemax:
    fig2b = plt.figure(figsize=(2*fw,1.25*fh))
    axmax = fig2b.gca()
    axmax.plot(t_range,upright['pax1']['md_list'],linestyle='dashed',lw=1,color=tab10(0),marker='o',mfc='None',mec=tab10(0),label='upright, ax1')
    axmax.plot(t_range,diag['pax1']['md_list'],linestyle='solid',lw=2,color=tab10(0),marker='o',mfc=tab10(0),mec=tab10(0),label='diag, ax1')
    axmax.plot(t_range,upright['pax2']['md_list'],linestyle='dashed',lw=1,color=tab10(2),marker='^',mfc='None',mec=tab10(2),label='upright,ax2')
    axmax.plot(t_range,diag['pax2']['md_list'],linestyle='solid',lw=2,color=tab10(2),marker='^',mfc=tab10(2),mec=tab10(2),label='diag,ax2')
    axmax.legend()
    axmax.set_xlabel('Time (hr)')
    axmax.set_ylabel(r'Max distance dye > $10^{-3}$')
    axmax.set_ylim([0,160])


if plotting:
    fig.subplots_adjust(right=0.98,left=0.09,bottom=0.05,top = 0.98,wspace=0.3,hspace=0.05)

plt.show(block=False)
plt.pause(0.1)
# plt.close()


