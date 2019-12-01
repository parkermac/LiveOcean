"""
Plot a basic salt budget for the segments.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
from datetime import datetime, timedelta

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zrfun

import tef_fun
import flux_fun

from time import time

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
# associated with lines like QQp[QQ<=0] = np.nan


# Get Ldir
Ldir = Lfun.Lstart('cas6', 'v3')

# select input/output location
indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

# load DataFrames of transport, volume, and salinity
v_df = pd.read_pickle(indir + 'volumes.p')
s_df = pd.read_pickle(indir + 'daily_segment_salinity.p')
# note s_df.index is daily datetimes

sect_df = tef_fun.get_sect_df()

seg_list = list(v_df.index)

vol = v_df.loc[seg_list,'volume m3'].sum()/1e9

# form a time series of the net salt in the system
net_s_df = pd.DataFrame(index=s_df.index, columns=['net_salt'])
ns = []
for dt in net_s_df.index:
     ns.append((s_df.loc[dt,seg_list]*v_df.loc[seg_list,'volume m3']).sum())
     # I had to do this as list items because if I added the net salt entries
     # to the DataFrame directly they ended up as objects of some sort, and then
     # later resample operations threw and error.
net_s_df.loc[:,'net_salt'] = ns
    
def get_fluxes(sect_name, in_sign=1):
    
    # form time series of net 2-layer transports into (+) and out of (-) the volume
    bulk = pickle.load(open(indir0 + 'bulk/' + sect_name + '.p', 'rb'))
    QQ = bulk['QQ']
    SS = bulk['SS']
    ot = bulk['ot']

    dt2 = []
    for tt in ot:
        dt2.append(Lfun.modtime_to_datetime(tt))

    # separate inflowing and outflowing transports
    QQp = QQ.copy()
    QQp[QQ<=0] = np.nan
    QQm = QQ.copy()
    QQm[QQ>=0] = np.nan

    # form two-layer versions of Q and S
    if in_sign == 1:
        Qin = np.nansum(QQp, axis=1)
        QSin = np.nansum(QQp*SS, axis=1)
        Qout = np.nansum(QQm, axis=1)
        QSout = np.nansum(QQm*SS, axis=1)
    elif in_sign == -1:
        Qin = -np.nansum(QQm, axis=1)
        QSin = -np.nansum(QQm*SS, axis=1)
        Qout = -np.nansum(QQp, axis=1)
        QSout = -np.nansum(QQp*SS, axis=1)
    Sin = QSin/Qin
    Sout = QSout/Qout

    tef_df = pd.DataFrame(index=dt2, columns=['Qin','Qout','QSin','QSout','Sin','Sout'])
    tef_df.loc[:,'Qin']=Qin
    tef_df.loc[:,'Qout']=Qout
    tef_df.loc[:,'QSin']=QSin
    tef_df.loc[:,'QSout']=QSout
    tef_df.loc[:,'Sin']=Sin
    tef_df.loc[:,'Sout']=Sout
    
    return tef_df
    
tef_df_jdf1 = get_fluxes('jdf1')
tef_df_sog5 = get_fluxes('sog5', in_sign=-1)

tef_df = pd.DataFrame(index=tef_df_jdf1.index, columns=['QSin','QSout','Snet'])
tef_df.loc[:,'QSin'] = tef_df_jdf1.loc[:,'QSin'] + tef_df_sog5.loc[:,'QSin']
tef_df.loc[:,'QSout'] = tef_df_jdf1.loc[:,'QSout'] + tef_df_sog5.loc[:,'QSout']

tef_mean_df = tef_df.resample('1M').mean()
# the above puts timestamps at the end of the month
# so here we set it to the middle of each month becasue it is more
# consistent with the averaging
tef_mean_df.index -= timedelta(days=15)
tef_mean_df.loc[:,'yd'] = tef_mean_df.index.dayofyear

# add monthly average net salt
for tdt in list(tef_mean_df.index):
    a = net_s_df[net_s_df.index.month==tdt.month].mean()
    tef_mean_df.loc[tdt,'Snet'] = float(a)
    
yd = np.array(tef_mean_df.loc[:,'yd'].values)
qs_in = np.array(tef_mean_df.loc[:,'QSin'].values)
qs_out = np.array(tef_mean_df.loc[:,'QSout'].values)
sn = np.array(tef_mean_df.loc[:,'Snet'].values)
dsdt = np.diff(sn) / (np.diff(yd)*86400)

# plotting
plt.close('all')
fig = plt.figure(figsize=(12,7))

ax = fig.add_subplot(111)
ax.plot(yd[1:], dsdt, '-r', yd[1:], qs_in[1:], '-b', yd[1:], -qs_out[1:], '-g',
    yd[1:], qs_in[1:]+qs_out[1:], '--r', linewidth=3)

ax.grid(True)

plt.show()
    

