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
import zfun

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

# load low passed segment volumne and net salt DataFrames
v_lp_df = pd.read_pickle(indir + 'daily_segment_volume.p')
sv_lp_df = pd.read_pickle(indir + 'daily_segment_net_salt.p')

# USER set which volume to consider
which_vol = 'Salish Sea'
if which_vol == 'Salish Sea':
    seg_list = list(v_lp_df.columns)
    sect_sign_dict = {'jdf1':1, 'sog5':-1}
elif which_vol == 'Puget Sound':
    seg_list = flux_fun.ssA + flux_fun.ssM + flux_fun.ssT + flux_fun.ssS + flux_fun.ssW + flux_fun.ssH
    sect_sign_dict = {'ai1':1, 'dp':1}
elif which_vol == 'Hood Canal':
    seg_list = flux_fun.ssH
    sect_sign_dict = {'hc1':1}

sv_lp_df = sv_lp_df[seg_list]

sect_df = tef_fun.get_sect_df()

sv_lp_df = sv_lp_df[seg_list]

river_list = []
for seg_name in seg_list:
    seg = flux_fun.segs[seg_name]
    river_list = river_list + seg['R']
riv_df = pd.read_pickle(Ldir['LOo'] + 'river/' + 'cas6_v3_2017.01.01_2017.12.31.p')
riv_df.index += timedelta(days=0.5)
riv_df = riv_df.loc[sv_lp_df.index, river_list]
    
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
    
tef_df_dict = {}
for sn in sect_sign_dict.keys():
    in_sign = sect_sign_dict[sn]
    tef_df_dict[sn] = get_fluxes(sn, in_sign=in_sign)

# volume budget
vol_df = pd.DataFrame(index=sv_lp_df.index, columns=['Qin','-Qout', 'Qr','dV_dt','dV_dt_check'])
vol_df.loc[:,['Qin','-Qout']] = 0
for sect_name in sect_sign_dict.keys():
    df = tef_df_dict[sect_name]
    vol_df.loc[:,'Qin'] = vol_df.loc[:,'Qin'] + df.loc[:,'Qin']
    vol_df.loc[:,'-Qout'] = vol_df.loc[:,'-Qout'] - df.loc[:,'Qout']
qr = riv_df.sum(axis=1)
vol_df.loc[:,'Qr'] = qr
v = v_lp_df[seg_list].sum(axis=1).values
vol_df.loc[1:-1, 'dV_dt'] = (v[2:] - v[:-2]) / (2*86400)
vol_df.loc[:,'dV_dt_check'] = vol_df.loc[:,'Qin'] - vol_df.loc[:,'-Qout'] + vol_df.loc[:,'Qr']
vol_err = (vol_df['dV_dt'] - vol_df['dV_dt_check']).std()
#vol_rel_err = vol_err/vol_df['Qin'].mean()

# alternate calculation: mean error / mean Qr
vol_rel_err = (vol_df['dV_dt'] - vol_df['dV_dt_check']).mean()/qr.mean()

# salt budget
salt_df = pd.DataFrame(index=sv_lp_df.index, columns=['QSin','-QSout','dSnet_dt','dSnet_dt_check'])
salt_df.loc[:,['QSin','-QSout']] = 0
for sect_name in sect_sign_dict.keys():
    df = tef_df_dict[sect_name]
    salt_df.loc[:,'QSin'] = salt_df.loc[:,'QSin'] + df.loc[:,'QSin']
    salt_df.loc[:,'-QSout'] = salt_df.loc[:,'-QSout'] - df.loc[:,'QSout']
sn = sv_lp_df[seg_list].sum(axis=1).values
salt_df.loc[1:-1, 'dSnet_dt'] = (sn[2:] - sn[:-2]) / (2*86400)
salt_df.loc[:,'dSnet_dt_check'] = salt_df.loc[:,'QSin'] - salt_df.loc[:,'-QSout']
salt_err = (salt_df['dSnet_dt'] - salt_df['dSnet_dt_check']).std()
salt_rel_err = salt_err/salt_df['QSin'].mean()

plt.close('all')
fig = plt.figure(figsize=(14,7))

ax = fig.add_subplot(121)
salt_df.plot(ax=ax, grid=True, title=which_vol + ' Salt Budget (g/kg m3/s)').legend(loc='center')
ax.text(.05,.9, 'RMSE / Mean QSin = %0.2f%%' % (salt_rel_err*100), transform=ax.transAxes, fontsize=14)

ax = fig.add_subplot(122)
vol_df.plot(ax=ax, grid=True, title='Volume Budget (m3/s)').legend(loc='center')
ax.text(.05,.9, 'Mean Error / Mean Qr = %0.2f%%' % (vol_rel_err*100), transform=ax.transAxes, fontsize=14)

plt.show()
