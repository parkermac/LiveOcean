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
    salt_df.loc[1:-1, 'dSnet_dt'] = ((sn[2:] - sn[:-2])/(2*86400))/1000
    salt_df.loc[:,'Smean'] = sm
    salt_df['Sin'] = salt_df['QSin'] / vol_df['Qin']
    salt_df['Sout'] = salt_df['-QSout'] / vol_df['-Qout']
    
    # combined
    combined_df = pd.DataFrame(index=sv_lp_df.index)
    combined_df['Smean'] = salt_df['Smean']
    combined_df['dSnet_dt'] = salt_df['dSnet_dt']
    combined_df['Sin'] = salt_df['Sin']
    combined_df['Sout'] = salt_df['Sout']
    combined_df['QSin'] = salt_df['QSin']/1000
    combined_df['-QSout'] = salt_df['-QSout']/1000
    combined_df['Qr'] = vol_df['Qr']/1000
    combined_df['Qin'] = vol_df['Qin']/1000
    combined_df['-Qout'] = vol_df['-Qout']/1000
    
    if yy == 0:
        all_years_df = combined_df.copy()
    else:
        all_years_df = all_years_df.append(combined_df)
        
    yy += 1
    

# add a few more columns to plot in a different way
all_years_df['Qe'] = (all_years_df['Qin'] + all_years_df['-Qout'])/2
all_years_df['Qrr'] = (all_years_df['-Qout'] - all_years_df['Qin'])
all_years_df['DS'] = all_years_df['QSin']/all_years_df['Qin'] - all_years_df['-QSout']/all_years_df['-Qout']
all_years_df['Sbar'] = (all_years_df['QSin']/all_years_df['Qin'] + all_years_df['-QSout']/all_years_df['-Qout'])/2
all_years_df['QeDS'] = all_years_df['Qe'] * all_years_df['DS']
all_years_df['-QrSbar'] = -all_years_df['Qrr'] * all_years_df['Sbar']

# Fix a problem where .resample() would drop the Sin column
# because it was somhow not numeric
for cn in all_years_df.columns:
    all_years_df[cn] = pd.to_numeric(all_years_df[cn])
# NOTE: I think .to_numeric() operates on Series, so we only do one column at a time.
all_years_df = all_years_df.resample('M', loffset='-15d').mean()

# plotting

plt.close('all')
fig = plt.figure(figsize=(14,8))

tx = .05
ty = .9
ty2 = .05
fs = 14
lw = 3
dt0 = datetime(2017,1,1)
dt1 = datetime(2020,1,1)

ax = fig.add_subplot(221)
all_years_df[['Smean', 'Sin','Sout']].plot(ax=ax, grid=True, color=['purple','r','orange'], linewidth=lw)
ax.legend(labels=['$S_{mean}$','$S_{in}$','$S_{out}$'], loc='lower right')
ax.text(tx, ty, '(a) ' + which_vol + ' Salinities $(g/kg)$', size=fs, transform=ax.transAxes)
ax.set_xticklabels([])
ax.set_xticklabels([], minor=True)
ax.set_xlim(dt0, dt1)

ax = fig.add_subplot(222)
all_years_df[['Qin', '-Qout']].plot(ax=ax, grid=True, color=['r','orange'], linewidth=lw)
ax.legend(labels=['$Q_{in}$','$-Q_{out}$'], loc='lower right')
ax.set_ylim(bottom=0)
ax.text(tx, ty2, '(b) Exchange Flow $(10^{3}\ m^{3}s^{-1})$', size=fs, transform=ax.transAxes)
ax.set_xticklabels([])
ax.set_xticklabels([], minor=True)
ax.set_xlim(dt0, dt1)

ax = fig.add_subplot(223)
all_years_df['Qr'].plot(ax=ax, grid=True, legend=False, color='c', linewidth=lw)
ax.set_ylim(bottom=0)
ax.text(tx, ty2, '(c) Net River Flow $(10^{3}\ m^{3}s^{-1})$', size=fs, transform=ax.transAxes)
ax.set_xlim(dt0, dt1)

ax = fig.add_subplot(224)
all_years_df[['dSnet_dt','QeDS', '-QrSbar']].plot(ax=ax, grid=True, color=['peru','b','g'], linewidth=lw)
ax.legend(labels=['$d S_{net} / dt$','$Q_e \Delta S$', '$-Q_R S_{bar}$'], loc='upper right')
ax.text(tx, ty, '(d) Salt Budget Terms $(g/kg\ 10^{3}\ m^{3}s^{-1})$', size=fs, transform=ax.transAxes)
ax.set_xlim(dt0, dt1)

fig.tight_layout()


plt.show()

