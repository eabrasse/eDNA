import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import numpy as np


plt.close('all')

home = '/Users/elizabethbrasseale/Projects/eDNA/'
decay_quant_fn = home+'data/All_Mod1_Molecular_Data - all_data.csv'
dfq = pd.read_csv(decay_quant_fn,sep=',',engine='python')
gq = dfq[dfq['campaign']=='Feb_2023']
df= gq[gq['task']=='decay']

df['datetime'] = pd.to_datetime(df.month.apply(str)+'/'+df.day.apply(str)+'/'+df.year.apply(str)+' '+df['filter_start_time '])
df['datetime'] = df['datetime'].dt.tz_localize(tz='America/Vancouver')
df['datetime'] = df['datetime'].dt.tz_convert('UTC')
df['timestamp'] = df['datetime'].view('int64').values//10**9
df['hrs_since_start'] = 1/3600*(df['timestamp']-df['timestamp'].min())
df=df.sort_values('hrs_since_start')


# C0 = df[df['hrs_since_start']==0]['dolphin_DNA_copy_uL'].mean()
C0 = df['dolphin_DNA_copy_uL'].max()
df['lnC_over_C0'] = df.apply(lambda row: np.log(row.dolphin_DNA_copy_uL/C0),axis=1)

fig, ax = plt.subplots()
ax.plot(df.hrs_since_start,df.lnC_over_C0,linestyle='none',marker='o',mfc='k',mec='k',markersize=3,zorder=100)
ax.set_xlabel('Hours since removal')
ax.set_ylabel('log(C/C0)')

df['hrs_since_start'] = df['hrs_since_start'].astype(float)
# X = sm.add_constant(df['hrs_since_start'].values)
# ols_model = sm.OLS(df['lnC_over_C0'].values, X)
ols_model = sm.OLS(df['lnC_over_C0'].values, df['hrs_since_start'].values)
est = ols_model.fit()
out = est.conf_int(alpha=0.05, cols=None)

# fig, ax = plt.subplots()
# df.plot(x='year',y='count',linestyle='None',marker='s', ax=ax)
# y_pred = est.predict(X)
y_pred = est.predict(df['hrs_since_start'].values)
x_pred = df.hrs_since_start.values
ax.plot(x_pred,y_pred,linestyle='solid',color='blue')
#
# pred = est.get_prediction(X).summary_frame()
pred = est.get_prediction(df['hrs_since_start'].values).summary_frame()
# ax.plot(x_pred,pred['mean_ci_lower'],linestyle='--',color='blue')
# ax.plot(x_pred,pred['mean_ci_upper'],linestyle='--',color='blue')
ax.fill_between(x_pred,pred['mean_ci_lower'],pred['mean_ci_upper'],color='gray',alpha=0.5,linewidth=0)
ax.grid()
#
# # Alternative way to plot
# def line(x,b=0,m=1):
#     return m*x+b
#
# ax.plot(x_pred,line(x_pred,est.params[0],est.params[1]),color='blue')

plt.show(block=False)
plt.pause(0.1)