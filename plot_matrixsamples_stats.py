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
import matplotlib.dates as mdates
from matplotlib.dates import DayLocator, HourLocator, MonthLocator,DateFormatter, drange

plt.close('all')
tab10 = plt.get_cmap('tab10',10)
tz = pytz.utc
pdt = pytz.timezone('America/Vancouver')

# load in data
home = '/Users/elizabethbrasseale/Projects/eDNA/'


ESP_fn = home+'data/MURI_Module1_datasheets - 03_ESP1_Feb2023.csv'
ESP_df = pd.read_csv(ESP_fn,sep=',',engine='python')
ESP_df['datetime'] = pd.to_datetime(ESP_df['ESP_date']+' '+ESP_df['ESP_filter_t0'])
ESP_df = ESP_df.dropna(subset=['ESP_date','PB_quantity'])
ESP_df['datetime'] = ESP_df['datetime'].dt.tz_localize(tz='America/Vancouver')
ESP_df['datetime'] = ESP_df['datetime'].dt.tz_convert('UTC')

datafn = home+'data/ESP_box_model_output.csv'
df = pd.read_csv(datafn,sep=',',engine='python')
df['Datetime']=pd.to_datetime(df['Datetime']).dt.tz_localize(tz='UTC')


## add matrix sample plotting code

#use Ryan's estimated eDNA quants
DNA_quant_fn = home+'data/projected_logDNA_concentrations.csv'
dfq = pd.read_csv(DNA_quant_fn,sep=',',engine='python')

# get metadata from master spreadsheet
master_fn = home+'data/All_Mod1_Molecular_Data - all_data.csv'
dfm = pd.read_csv(master_fn,sep=',',engine='python')

#create a unique identifier per sample including lab ID, biorep, and technrep
for df0 in dfq,dfm:
    df0['data_ID'] = df0.apply(lambda row: row.lab_ID + str(row.bio_rep) + str(row.tech_rep),axis=1)

# just look at matrix samples
dfm = dfm[dfm['task']=='matrix']
dfm = dfm[~dfm['reference'].str.startswith('SR_later_PCI_')]
# dfm.loc[dfm['collection_time']=='18:30','collection_time']='15:00' # fix inaccurate value in spreadsheet
dfm = dfm[dfm['collection_time']!='18:30']

#DNA quants to master dataframe
dfm['mean_est']=dfq[dfq.data_ID.isin(dfm.data_ID)].mean_est.values
dfm['p0.25']=dfm['mean_est']-dfq[dfq.data_ID.isin(dfm.data_ID)]['p0.25'].values
dfm['p0.75']=dfq[dfq.data_ID.isin(dfm.data_ID)]['p0.75'].values-dfm['mean_est']

# dfm['dist_from_pen'] = dfm.apply(lambda row: efun.ll2dist(row.lon,row.lat,lon0,lat0),axis=1) 
dfm['hour'] = dfm.apply(lambda row: row.collection_time[:2],axis=1)
dfm['minute'] = dfm.apply(lambda row: row.collection_time[3:],axis=1)
# dfm['datetime'] = pd.to_datetime('{}-{}-{} {}:00'.format(str(dfm.year),str(dfm.month),str(dfm.day),str(dfm.collection_time)))
dfm['datetime'] = pd.to_datetime({'year':dfm.year,'month':dfm.month,'day':dfm.day,'hour':dfm.hour,'minute':dfm.minute})
dfm['timestamp'] = dfm['datetime'].astype(int)
dfm['seconds_since_t0'] = dfm.timestamp-dfm.timestamp.values[0]
dfm['seconds_since_t0'] = dfm.seconds_since_t0/10**9

fw,fh = efun.gen_plot_props()


extraction_list = dfm['extraction'].unique()
preservation_list = dfm['preservation'].unique()

for extraction in extraction_list:
    for preservation in preservation_list:
        gg = dfm[(dfm.extraction==extraction)&(dfm.preservation==preservation)]

        if np.shape(gg)[0]>0:

            fig = plt.figure(figsize=(2*fw,fh))
            ax = fig.gca()

            ax.plot(gg.datetime,np.log(gg.dolphin_DNA_copy_uL),linestyle='none',marker='o',mfc='None',mec=tab10(0),markersize=8)
            ax.set_ylabel(r'log DNA conc (copies $\mu\mathrm{L}^{-1}$)')
            ax.xaxis.set_major_locator(mdates.MinuteLocator(range(0,60,10)))

            ax.set_xlabel('time')
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            
            techrepmeans = []
            techrepvars = []
            techrepcvs = []
            for bio_rep in gg.bio_rep.unique():
                ggg= gg[gg.bio_rep==bio_rep]
                trm = np.mean(np.log(ggg.dolphin_DNA_copy_uL))
                trs = np.std(np.log(ggg.dolphin_DNA_copy_uL))
                trc = trs/trm
                techrepmeans.append(trm)
                techrepvars.append(trs**2)
                techrepcvs.append(trc)
                
            mean = np.mean(techrepmeans)
            std = np.std(techrepmeans)
            cv = std/mean
            var = std**2
            
            techvar = np.mean(techrepvars)
            
            ax.text(0.2,0.9,f'{extraction:}, {preservation:}\nmean = {mean:2.2E}\nvar = {var:2.2E}\nCV = {cv:.2E}\ntech var = {techvar:.2E}', transform=ax.transAxes,va='top')
            # print(f'{extraction:}, {preservation:}')

            # fig.show()
            fig.subplots_adjust(left=0.2,right=0.98,bottom=0.2)
            # fig.pause(0.1)
            outfn = home+ f'plots/matrix stats/LOGNORMAL_ESP_and_matrix_samples_{extraction:}_{preservation:}.png'
            plt.savefig(outfn,format='png',dpi=600)

# extraction = 'QBTK'
# preservation = 'RNAlater'
# gg = dfm[(dfm.extraction==extraction)&(dfm.preservation==preservation)]
# mean = np.mean(gg.dolphin_DNA_copy_uL)
# std = np.std(gg.dolphin_DNA_copy_uL)