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
which_vol = 'Puget Sound'

if which_vol == 'Salish Sea':
    seg_list = list(v_lp_df.columns)
    seg_list = seg_list[1:]
    sect_sign_dict = {'jdf2':1, 'sog5':-1}
elif which_vol == 'Puget Sound':
    seg_list = flux_fun.ssA + flux_fun.ssM + flux_fun.ssT + flux_fun.ssS + flux_fun.ssW + flux_fun.ssH
    seg_list = seg_list[1:]
    sect_sign_dict = {'ai2':1, 'dp':1}
elif which_vol == 'Hood Canal':
    seg_list = flux_fun.ssH
    seg_list = seg_list[1:]
    sect_sign_dict = {'hc2':1}
elif which_vol == 'South Sound':
    seg_list = flux_fun.ssT + flux_fun.ssS
    seg_list = seg_list[1:]
    sect_sign_dict = {'tn2':1}

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
    fnet = bulk['fnet_lp'] # net tidal energy flux
    tef_df = pd.DataFrame(index=dt2, columns=['Qin','Qout','QSin','QSout','Sin','Sout','fnet'])
    tef_df.loc[:,'Qin']=Qin
    tef_df.loc[:,'Qout']=Qout
    tef_df.loc[:,'QSin']=QSin
    tef_df.loc[:,'QSout']=QSout
    tef_df.loc[:,'Sin']=Sin
    tef_df.loc[:,'Sout']=Sout
    tef_df.loc[:,'fnet'] = fnet * in_sign

    return tef_df
    
tef_df_dict = {}
for sn in sect_sign_dict.keys():
    in_sign = sect_sign_dict[sn]
    tef_df_dict[sn] = get_fluxes(sn, in_sign=in_sign)
    
# hack to get dsdx
if which_vol == 'Salish Sea':
    df1 = get_fluxes('jdf1')
    df3 = get_fluxes('jdf3')
elif which_vol == 'Puget Sound':
    df1 = get_fluxes('ai1')
    df3 = get_fluxes('ai3')
elif which_vol == 'South Sound':
    df1 = get_fluxes('tn1')
    df3 = get_fluxes('tn3')
elif which_vol == 'Hood Canal':
    df1 = get_fluxes('hc1')
    df3 = get_fluxes('hc3')
dsdx = (df1['Sin']+df1['Sout'])/2-(df3['Sin']+df3['Sout'])/2
dissip = df1['fnet']-df3['fnet']

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
salt_df = pd.DataFrame(index=sv_lp_df.index, columns=['QSin','-QSout','dSnet_dt','dSnet_dt_check','fnet'])
salt_df.loc[:,['QSin','-QSout','fnet']] = 0
for sect_name in sect_sign_dict.keys():
    df = tef_df_dict[sect_name]
    salt_df.loc[:,'QSin'] = salt_df.loc[:,'QSin'] + df.loc[:,'QSin']
    salt_df.loc[:,'-QSout'] = salt_df.loc[:,'-QSout'] - df.loc[:,'QSout']
    salt_df.loc[:,'fnet'] = salt_df.loc[:,'fnet'] + df.loc[:,'fnet']
sn = sv_lp_df[seg_list].sum(axis=1).values
salt_df.loc[1:-1, 'dSnet_dt'] = (sn[2:] - sn[:-2]) / (2*86400)
salt_df.loc[:,'dSnet_dt_check'] = salt_df.loc[:,'QSin'] - salt_df.loc[:,'-QSout']
salt_err = (salt_df['dSnet_dt'] - salt_df['dSnet_dt_check']).std()
salt_rel_err = salt_err/salt_df['QSin'].mean()
# add a few more columns to plot in a different way
salt_df['Qe'] = (vol_df['Qin'] + vol_df['-Qout'])/2
salt_df['Qr'] = (vol_df['-Qout'] - vol_df['Qin'])
salt_df['DS'] = salt_df['QSin']/vol_df['Qin'] - salt_df['-QSout']/vol_df['-Qout']
salt_df['Sbar'] = (salt_df['QSin']/vol_df['Qin'] + salt_df['-QSout']/vol_df['-Qout'])/2
salt_df['QeDS'] = salt_df['Qe'] * salt_df['DS']
salt_df['-QrSbar'] = -salt_df['Qr'] * salt_df['Sbar']

sf_df = salt_df[['QeDS','-QrSbar','dSnet_dt']]
sf_df['Error'] = sf_df['dSnet_dt'] - sf_df['QeDS'] - sf_df['-QrSbar']
sf_df['Forcing'] = (dsdx**3 / dissip)
sf_df['Dissipation'] = dissip
sf_df['dSdx'] = dsdx

plt.close('all')
fig = plt.figure(figsize=(20,10))

#ax = plt.subplot2grid((1,3), (0,0), colspan=2)
ax = fig.add_subplot(211)
sf_df[['QeDS','-QrSbar','dSnet_dt','Error']].plot(ax=ax, grid=True,
    title=which_vol + ' Salt Budget (g/kg m3/s)', legend=False)
legh = ax.legend(labels=['$Q_e \Delta S$', '$-Q_R S_{bar}$', '$d S_{net} / dt$','Error'])
ii = 4
for frc in ['Forcing', 'Dissipation', 'dSdx']:
    ax = fig.add_subplot(2,3,ii)
    sf_df.plot(x=frc, y='QeDS', ax=ax, grid=True, style='og', markersize=2)
    ax.set_xlabel('')
    ax.text(.05,.9,frc, transform=ax.transAxes)
    ii += 1

plt.show()
