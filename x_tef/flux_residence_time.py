"""
Designed to calculate residence time from any of the "ic" experiments.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
import zfun

import tef_fun
import flux_fun

# required command line arguments, can be input in any order
parser = argparse.ArgumentParser()
parser.add_argument('-src', '--source', nargs='?', type=str, default='ic_hood_canal_inner')
args = parser.parse_args()
source = args.source

print(source)
if 'ic_' not in source:
    print(' -- Error: neet to run with an ic_ source --')
    sys.exit()

# Input directory
indir0 = Ldir['LOo'] + 'tef/'
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/flux/'

# hacky way of getting the year, assumes "item" is of the form:
# 'cas6_v3_lo8b_2017.01.01_2017.12.31'
year_str = item.split('_')[-1].split('.')[0]
year = int(year_str)

print(item)

# load a Series of the volumes of each segment, created by flux_get_vol.py
v_df = pd.read_pickle(indir + 'volumes.p')
V = flux_fun.get_V(v_df)

for season in flux_fun.season_list:
    
    seg2_list = flux_fun.ic_seg2_dict[source]
        
    infile = 'aa_' + source + '_' + season + '.p'

    aa = pd.read_pickle(indir + infile)
    
    this_aa = aa.loc[:,seg2_list]
    this_V = V[seg2_list]
    net_V = this_V.sum()

    this_net_aa = this_aa.copy()

    for sn in this_V.index:
        VV = this_V[sn]
        this_net_aa.loc[:,sn] = this_net_aa.loc[:,sn] * VV
    
    # make a series of mean concentration in the volume
    mean_c = pd.Series(index=this_aa.index)
    mean_c = this_net_aa.sum(axis=1) / net_V

    # find e-folding time
    td = mean_c.index.values
    mc = mean_c.values

    ind_ef = np.argwhere(mc < 1/np.e)[0]
    
    print(' %10s: Tres = %0.2f days' % (season, td[ind_ef]))

    if False:
        fig = plt.figure(figsize=(13,8))
        ax = fig.add_subplot(111)
        fs = 16
        lw = 3
        mean_c.plot(ax=ax, grid=True, linewidth=lw)
        ttext = infile.replace('ic_','').replace('aa_','').replace('.p','').replace('_',' ').title()
        ax.set_title(season.title() + ' -- IC: ' + ttext, fontsize=fs)
        ax.set_xlabel('Time (days)', fontsize=fs)
        ax.set_ylabel('Mean Concentration in Volume', fontsize=fs)
        ax.set_ylim(0,1)
        ax.tick_params(labelsize=fs-2) 
        ax.plot(td[ind_ef], mc[ind_ef], '*r', markersize=fs+3)
        ax.text(.95,.85,'%s e-Folding Time = %d Days' % (year_str, td[ind_ef]), fontsize=fs+2,
            transform=ax.transAxes, horizontalalignment='right', color='r', fontweight='bold')
        plt.show()
        outdir = indir0 + item + '/misc_figs/'
        #fig.savefig(outdir + 'HC_residence_time_plot_' + season + '.png')

