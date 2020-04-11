"""
Calculates an exchange-flow oriented salt budget from TEF terms, and
explores dynamical scaling:

does Qe behave as expected relative to dSbar_dx and K?

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
from importlib import reload
reload(flux_fun)

from time import time

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', type=str, default='cas6')
parser.add_argument('-t', '--tag', type=str, default='v3')
parser.add_argument('-x', '--ex_name', type=str, default='lo8b')
parser.add_argument('-y', '--year', type=int, default=2017)
parser.add_argument('-v', '--volume', type=str, default='Puget Sound')
args = parser.parse_args()
which_vol = args.volume
year_str = str(args.year)

# Get Ldir
Ldir = Lfun.Lstart(args.gridname, args.tag)
gtagex = args.gridname + '_' + args.tag + '_' + args.ex_name

# select input/output location
run_name = gtagex+'_'+year_str+'.01.01_'+year_str+'.12.31'
indir00 = Ldir['LOo'] + 'tef/'
indir0 = indir00 + run_name + '/'
indir = indir0 + 'flux/'

# load low passed segment volume and net salt DataFrames
v_lp_df = pd.read_pickle(indir + 'daily_segment_volume.p')
sv_lp_df = pd.read_pickle(indir + 'daily_segment_net_salt.p')

# get volumes
voldir = indir00 + 'volumes_' + Ldir['gridname'] + '/'
v_df = pd.read_pickle(voldir + 'volumes.p')

# get section definitions
sect_df = tef_fun.get_sect_df()

if which_vol == 'Salish Sea':
    seg_list = list(v_lp_df.columns)
    sect_sign_dict = {'jdf1':1, 'sog5':-1}
elif which_vol == 'Puget Sound':
    seg_list = (flux_fun.ssA + flux_fun.ssM + flux_fun.ssT
        + flux_fun.ssS + flux_fun.ssW + flux_fun.ssH)
    sect_sign_dict = {'ai1':1, 'dp':1}
elif which_vol == 'Hood Canal':
    seg_list = flux_fun.ssH
    sect_sign_dict = {'hc1':1}
elif which_vol == 'South Sound':
    seg_list = flux_fun.ssT + flux_fun.ssS
    sect_sign_dict = {'tn1':1}

v_lp_df = v_lp_df[seg_list]
sv_lp_df = sv_lp_df[seg_list]

river_list = []
for seg_name in seg_list:
    seg = flux_fun.segs[seg_name]
    river_list = river_list + seg['R']
    riv_df = pd.read_pickle(Ldir['LOo'] + 'river/'
        + Ldir['gtag'] + '_'+year_str+'.01.01_'+year_str+'.12.31.p')
riv_df.index += timedelta(days=0.5)
riv_df = riv_df[river_list]
        
tef_df_dict = {}
for sn in sect_sign_dict.keys():
    in_sign = sect_sign_dict[sn]
    tef_df_dict[sn] = flux_fun.get_fluxes(indir0, sn, in_sign=in_sign)
    
vol_df, salt_df, vol_rel_err, salt_rel_err, salt_rel_err_qe = flux_fun.get_budgets(
    sv_lp_df, v_lp_df, riv_df, tef_df_dict, seg_list)

plt.close('all')
fig = plt.figure(figsize=(15,10))

ax = fig.add_subplot(111)
salt_df[['dSnet_dt','-QrSbar','QeDS','Error']].plot(ax=ax, grid=True, legend=False)
ax.set_title('%s Salt Budget $(g\ kg^{-1}\ m^{3}s^{-1})$' % (which_vol))
ax.set_xlim(salt_df.index[0], salt_df.index[-1])
legh = ax.legend(labels=['$d S_{net} / dt$', '$Q_{net} S_{bar}$', '$Q_e \Delta S$', 'Error'])
ax.text(.05,.9, 'Mean Error / Mean QeDS = %0.2f%%' % (salt_rel_err_qe*100), transform=ax.transAxes, fontsize=14)

plt.show()
