"""
The goal here is to see if there is any relationship between the
average salinity of a specified volume and a filtered version of
the forcing by rivers or QSin.

Designed to run over three years, so that we capture the effect
of the increasing salinity from 2017 to 2019.

"""

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import tef_fun
import flux_fun

import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
from datetime import datetime, timedelta
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', type=str, default='cas6')
parser.add_argument('-t', '--tag', type=str, default='v3')
parser.add_argument('-x', '--ex_name', type=str, default='lo8b')
parser.add_argument('-v', '--volume', type=str, default='Puget Sound')
args = parser.parse_args()
which_vol = args.volume

# Get Ldir
Ldir = Lfun.Lstart(args.gridname, args.tag)
gtagex = args.gridname + '_' + args.tag + '_' + args.ex_name

def get_fluxes(sect_name, in_sign=1):
    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages
    # associated with lines like QQp[QQ<=0] = np.nan
    
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

yy = 0
for year in [2017, 2018, 2019]:
    year_str = str(year)
    
    # select input/output location
    run_name = gtagex+'_'+year_str+'.01.01_'+year_str+'.12.31'
    indir00 = Ldir['LOo'] + 'tef/'
    indir0 = indir00 + run_name + '/'

    # outdir = indir00 + 'salt_budget_plots/'
    # Lfun.make_dir(outdir)
    # outname = outdir + 'salt_budget_' + year_str + '_' + which_vol.replace(' ','_') + '.png'

    # load low passed segment volume and net salt DataFrames
    v_lp_df = pd.read_pickle(indir0 + 'flux/daily_segment_volume.p')
    sv_lp_df = pd.read_pickle(indir0 + 'flux/daily_segment_net_salt.p')

    # info specific to each volume
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
    v_lp_df = v_lp_df[seg_list]
    
    river_list = []
    for seg_name in seg_list:
        seg = flux_fun.segs[seg_name]
        river_list = river_list + seg['R']
    riv_df = pd.read_pickle(Ldir['LOo'] + 'river/' + Ldir['gtag'] + '_'+year_str+'.01.01_'+year_str+'.12.31.p')
    riv_df.index += timedelta(days=0.5)
    riv_df = riv_df.loc[sv_lp_df.index, river_list]
    
    tef_df_dict = {}
    for sn in sect_sign_dict.keys():
        in_sign = sect_sign_dict[sn]
        tef_df_dict[sn] = get_fluxes(sn, in_sign=in_sign)

    # volume
    vol_df = pd.DataFrame(index=sv_lp_df.index, columns=['Qin','-Qout', 'Qr','V'])
    vol_df.loc[:,['Qin','-Qout']] = 0
    for sect_name in sect_sign_dict.keys():
        df = tef_df_dict[sect_name]
        vol_df.loc[:,'Qin'] = vol_df.loc[:,'Qin'] + df.loc[:,'Qin']
        vol_df.loc[:,'-Qout'] = vol_df.loc[:,'-Qout'] - df.loc[:,'Qout']
    qr = riv_df.sum(axis=1)
    vol_df.loc[:,'Qr'] = qr
    v = v_lp_df.sum(axis=1).values
    vol_df.loc[:,'V'] = v

    # salt
    salt_df = pd.DataFrame(index=sv_lp_df.index, columns=['QSin','-QSout','Snet','Smean', 'Sin'])
    salt_df.loc[:,['QSin','-QSout']] = 0
    for sect_name in sect_sign_dict.keys():
        df = tef_df_dict[sect_name]
        salt_df.loc[:,'QSin'] = salt_df.loc[:,'QSin'] + df.loc[:,'QSin']
        salt_df.loc[:,'-QSout'] = salt_df.loc[:,'-QSout'] - df.loc[:,'QSout']
    sn = sv_lp_df.sum(axis=1).values
    sm = sn/v
    salt_df.loc[:, 'Snet'] = sn
    salt_df.loc[:,'Smean'] = sm
    salt_df['Sin'] = salt_df['QSin'] / vol_df['Qin']
    
    # combined
    combined_df = pd.DataFrame(index=sv_lp_df.index, columns=['Qr','Smean','Sin'])
    combined_df['Smean'] = salt_df['Smean']
    combined_df['Sin'] = salt_df['Sin']
    combined_df['Qr'] = vol_df['Qr']
    
    if yy == 0:
        all_years_df0 = combined_df.copy()
        print(all_years_df0.columns)
    else:
        all_years_df0 = all_years_df0.append(combined_df)
        
    yy += 1
    
# Fix a problem where .resample() would drop the Sin column
# because it was somhow not numeric
all_years_df0['Sin'] = pd.to_numeric(all_years_df0['Sin'])
# NOTE: I think .to_numeric() operates on Series, so we only do one
# column at a time.
all_years_df = all_years_df0.resample('M').mean()

qrr = all_years_df['Qr'].to_numpy()
qrr_lagged = qrr.copy() * np.nan
lag = 2
qrr_lagged[lag:] = qrr[:-lag]
all_years_df.loc[:,'Qr_lagged'] = qrr_lagged

plt.close('all')
fig = plt.figure(figsize=(14,7))

ax = fig.add_subplot(221)
all_years_df['Smean'].plot(ax=ax, legend=False, grid=True)
ax.set_ylabel('Smean')

ax = fig.add_subplot(223)
all_years_df[['Qr', 'Qr_lagged']].plot(ax=ax, grid=True)
ax.set_ylabel('Qr')

ax = fig.add_subplot(222)
all_years_df.plot(x='Qr', y = 'Smean', style='ob', ax=ax, label='regular')
all_years_df.plot(x='Qr_lagged', y = 'Smean', style='or', ax=ax, label='lagged by ' + str(lag) + ' months')
ax.set_xlabel('Qr')
ax.set_ylabel('Smean')
ax.set_title(which_vol)

ax = fig.add_subplot(224)
all_years_df.plot(x='Sin', y = 'Smean', style='*g', ax=ax)

plt.show()

