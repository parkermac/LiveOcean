"""
Do salt and volume budgets for user-specified volumes.
This is a fairly complicated task because it involves coordinating
the net salt and volume time series from flux_[get,lowpass]_s.py,
the TEF transports through any open sections of the volume,
and all the rivers flowing into the volume.

Run with a command like:
run flux_salt_budget
and then it runs for all three years and all three volumes (Salish, PS, HC).

Result: typically the errors are acceptably small, like 1%
depending on how you define "error".

"""

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import tef_fun
import flux_fun
from importlib import reload
reload(flux_fun)

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
#parser.add_argument('-v', '--volume', type=str, default='Puget Sound')
args = parser.parse_args()
#which_vol = args.volume

# Get Ldir
Ldir = Lfun.Lstart(args.gridname, args.tag)
gtagex = args.gridname + '_' + args.tag + '_' + args.ex_name

# vol_list = ['Salish Sea', 'Puget Sound', 'Hood Canal']
# year_list = [2017, 2018, 2019]
vol_list = ['Puget Sound']
year_list = [2017]

err_df_vol = pd.DataFrame(index=year_list, columns=vol_list)
err_df_salt = pd.DataFrame(index=year_list, columns=vol_list)

plt.close('all')
for which_vol in vol_list:

    for year in year_list:
        year_str = str(year)
    
        # select input/output location
        run_name = gtagex+'_'+year_str+'.01.01_'+year_str+'.12.31'
        indir00 = Ldir['LOo'] + 'tef/'
        indir0 = indir00 + run_name + '/'

        outdir = indir00 + 'salt_budget_plots/'
        Lfun.make_dir(outdir)
        outname = outdir + 'salt_budget_' + year_str + '_' + which_vol.replace(' ','_') + '.png'

        # load low passed segment volume and net salt DataFrames
        v_lp_df = pd.read_pickle(indir0 + 'flux/daily_segment_volume.p')
        sv_lp_df = pd.read_pickle(indir0 + 'flux/daily_segment_net_salt.p')

        # info specific to each volume
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

        v_lp_df = v_lp_df[seg_list]
        sv_lp_df = sv_lp_df[seg_list]
    
        sect_df = tef_fun.get_sect_df()

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
            tef_df_dict[sn] =flux_fun.get_fluxes(indir0, sn, in_sign=in_sign)
        
        vol_df, salt_df, vol_rel_err, salt_rel_err, salt_rel_err_qe = flux_fun.get_budgets(
            sv_lp_df, v_lp_df, riv_df, tef_df_dict, seg_list)

        fig = plt.figure(figsize=(14,7))

        ax = fig.add_subplot(121)
        salt_df[['dSnet_dt','QSin','-QSout']].plot(ax=ax, grid=True)
        ax.legend(labels=[r'$Storage^{ADJ}$',r'$Q_{in}\Delta S$',r'$-Q_{R}S_{out}$',r'$Error$']loc='upper right')
        ax.set_title(year_str + ' ' + which_vol + ' Salt Budget (g/kg m3/s)')
        ax.text(.05,.9, 'Mean Error / Mean QSin = %0.2f%%' % (salt_rel_err*100), transform=ax.transAxes, fontsize=14)

        ax = fig.add_subplot(122)
        vol_df[['dV_dt','Qin','-Qout', 'Qr']].plot(ax=ax, grid=True).legend(loc='upper right')
        ax.set_title('Volume Budget (m3/s)')
        ax.text(.05,.9, 'Mean Error / Mean Qr = %0.2f%%' % (vol_rel_err*100), transform=ax.transAxes, fontsize=14)

        plt.savefig(outname)
        
        # save error statistics (all %)
        err_df_vol.loc[year, which_vol] = vol_rel_err*100
        err_df_salt.loc[year, which_vol] = salt_rel_err*100
        
plt.show()

