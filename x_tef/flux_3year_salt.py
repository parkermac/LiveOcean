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

indir00 = Ldir['LOo'] + 'tef/'
outdir = indir00 + 'misc_figs_cas6/'

qr_dict = {} # save annual mean river flow
yy = 0
for year in [2017, 2018, 2019]:
    year_str = str(year)
    
    # select input/output location
    run_name = gtagex+'_'+year_str+'.01.01_'+year_str+'.12.31'
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
    riv_df = pd.read_pickle(Ldir['LOo'] + 'river/'
        + Ldir['gtag'] + '_'+year_str+'.01.01_'+year_str+'.12.31.p')
    riv_df.index += timedelta(days=0.5)
    riv_df = riv_df.loc[sv_lp_df.index, river_list]
    
    tef_df_dict = {}
    for sn in sect_sign_dict.keys():
        in_sign = sect_sign_dict[sn]
        tef_df_dict[sn] = flux_fun.get_fluxes(indir0, sn, in_sign=in_sign)

    vol_df, salt_df, vol_rel_err, salt_rel_err, salt_rel_err_qe = flux_fun.get_budgets(
        sv_lp_df, v_lp_df, riv_df, tef_df_dict, seg_list)
        
    qr_dict[year] = vol_df['Qr'].mean()
        
    if yy == 0:
        vol_df_all = vol_df.copy()
        salt_df_all = salt_df.copy()
    else:
        vol_df_all = vol_df_all.append(vol_df)
        salt_df_all = salt_df_all.append(salt_df)
        
    yy += 1
    
# Fix a problem where .resample() would drop the Sin column
# because it was somhow not numeric
for cn in vol_df_all.columns:
    vol_df_all[cn] = pd.to_numeric(vol_df_all[cn])
for cn in salt_df_all.columns:
    salt_df_all[cn] = pd.to_numeric(salt_df_all[cn])
# NOTE: I think .to_numeric() operates on Series, so we only do one column at a time.
vol_df_all = vol_df_all.resample('M', loffset='-15d').mean()
salt_df_all = salt_df_all.resample('M', loffset='-15d').mean()

# rescale to be 1000 m3/s
for vn in ['QSin', '-QSout', 'Ftide', 'dSnet_dt', 'Error', 'Qe', 'Qnet', 'QeDS', '-QrSbar']:
    salt_df_all[vn] = salt_df_all[vn]/1000
for vn in ['Qin', '-Qout', 'Qtide', 'Qr', 'V', 'dV_dt', 'Error']:
    vol_df_all[vn] = vol_df_all[vn]/1000

# plotting


tx = .05
ty = .9
ty2 = .05
fs = 16
lw = 3
dt0 = datetime(2017,1,1)
dt1 = datetime(2020,1,1)

plt.rc('font', size=fs)
plt.close('all')
plt.close('all')
fig = plt.figure(figsize=(18,10))


ax = fig.add_subplot(221)
salt_df_all[['Smean', 'Sin','Sout']].plot(ax=ax, grid=False, color=['goldenrod','r','b'], linewidth=lw)
ax.legend(labels=['$S_{mean}$','$S_{in}$','$S_{out}$'], loc='lower right')
ax.text(tx, ty, '(a) ' + which_vol + ' Salinities $[g \ kg^{-1}]$', size=fs, transform=ax.transAxes,
    bbox=dict(facecolor='w', edgecolor='None',alpha=.5), weight='bold')
ax.set_xticklabels([])
ax.set_xticklabels([], minor=True)
ax.set_xlim(dt0, dt1)
ax.set_ylim(29,33)
ax.set_yticks([30,31,32])
ax.set_yticklabels([30,31,32])
#
ax.vlines([datetime(2018,1,1),datetime(2019,1,1)],29,33, alpha=.5)
ax.set_xticks([datetime(2017,1,1),datetime(2017,7,1),datetime(2018,1,1),
    datetime(2018,7,1),datetime(2019,1,1),datetime(2019,7,1),datetime(2019,12,31)])

ax = fig.add_subplot(222)
vol_df_all[['Qin', '-Qout']].plot(ax=ax, grid=False, color=['r','b'], linewidth=lw)
ax.legend(labels=['$Q_{in}$','$-Q_{out}$'], loc='lower right')
ax.set_ylim(bottom=0)
ax.text(tx, ty2, '(b) Exchange Flow $[1000\ m^{3}s^{-1}]$', size=fs, transform=ax.transAxes,
    bbox=dict(facecolor='w', edgecolor='None',alpha=.5), weight='bold')
ax.set_xticklabels([])
ax.set_xticklabels([], minor=True)
ax.set_xlim(dt0, dt1)
ax.set_ylim(30,50)
ax.set_yticks([30,40,50])
ax.set_yticklabels([30,40,50])
#
ax.vlines([datetime(2018,1,1),datetime(2019,1,1)],30,50, alpha=.5)
ax.set_xticks([datetime(2017,1,1),datetime(2017,7,1),datetime(2018,1,1),
    datetime(2018,7,1),datetime(2019,1,1),datetime(2019,7,1),datetime(2019,12,31)])

ax = fig.add_subplot(223)
vol_df_all['Qr'].plot(ax=ax, grid=False, legend=False, color='c', linewidth=lw)
ax.set_ylim(bottom=0)
ax.text(tx, ty2, '(c) Net River Flow $[1000 \ m^{3}s^{-1}]$', size=fs, transform=ax.transAxes,
    bbox=dict(facecolor='w', edgecolor='None',alpha=.5), weight='bold')
ax.set_xlim(dt0, dt1)
ax.set_ylim(0,3)
ax.set_yticks([0,1,2])
ax.set_yticklabels([0,1,2])
#
ax.vlines([datetime(2018,1,1),datetime(2019,1,1)],0,15, alpha=.5)
ax.set_xticks([datetime(2017,1,1),datetime(2017,7,1),datetime(2018,1,1),
    datetime(2018,7,1),datetime(2019,1,1),datetime(2019,7,1),datetime(2019,12,31)])
ax.set_xticklabels(['','2017','','2018','','2019',''], rotation=0,
    fontdict={'horizontalalignment':'center'})
#add annual means
for year in [2017,2018,2019]:
    ax.text(datetime(year,7,1), 2.5,
        'Mean =\n' + str(int(qr_dict[year])) + ' $[m^{3}s^{-1}]$',
        ha='center', va='center',
        color = 'c', weight='bold')

ax = fig.add_subplot(224)
salt_df_all[['dSnet_dt','QeDS', '-QrSbar']].plot(ax=ax, grid=False,
    color=['sandybrown','darkorchid','cornflowerblue'],
    linewidth=lw)
ax.legend(labels=['$d S_{net} / dt$','$Q_e \Delta S$', '$-Q_R S_{bar}$'], loc='upper right')
ax.text(tx, ty, '(d) Salt Budget Terms $[1000 \ g \ kg^{-1} \ m^{3}s^{-1}]$',
    size=fs, transform=ax.transAxes,
    bbox=dict(facecolor='w', edgecolor='None',alpha=.5), weight='bold')
ax.set_xlim(dt0, dt1)
ax.set_ylim(-100,100)
#
ax.vlines([datetime(2018,1,1),datetime(2019,1,1)],-100,100, alpha=.5)
ax.hlines(0,dt0,dt1)
ax.set_xticks([datetime(2017,1,1),datetime(2017,7,1),datetime(2018,1,1),
    datetime(2018,7,1),datetime(2019,1,1),datetime(2019,7,1),datetime(2019,12,31)])
ax.set_xticklabels(['','2017','','2018','','2019',''], rotation=0,
    fontdict={'horizontalalignment':'center'})

fig.tight_layout()
fig.savefig(outdir + 'PS_3year.png')
plt.show()
plt.rcdefaults()


