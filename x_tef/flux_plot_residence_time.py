"""
Plot the results of the flux age engine.  Designed to calculate
residence time from any of the "ic" experiments.

NOTE: Currently hardwired to only do the Hood Canal experiment.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
import zfun

import tef_fun
import flux_fun

# Input directory
indir0 = Ldir['LOo'] + 'tef/'
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/flux/'

outdir = indir0 + item + '/misc_figs/'

plt.close('all')
for season in ['full']:#flux_fun.season_list:
    
    # load the DataFrame of time series results of flux_engine.py
    # infile = Lfun.choose_item(indir, tag='aa_', itext='Choose flux engine output file:')
    # if 'ic_' not in infile:
    #     print('Please choose a file from an "ic_" experiment.')
    #     sys.exit()

    infile = Lfun.choose_item(indir, tag='aa_ic')
    
    #infile = 'aa_ic_hood_canal_' + season + '.p'

    aa = pd.read_pickle(indir + infile)

    # load a Series of the volumes of each segment, created by flux_get_vol.py
    v_df = pd.read_pickle(indir + 'volumes.p')

    # create Series of two-layer volumes
    # Since this HAS to be the same as what the flux_engine used, it really should
    # be put in a function somewhere else.
    V = pd.Series(index=v_df.index)
    for seg_name in v_df.index:
        V[seg_name+'_s'] = 0.8 * v_df.loc[seg_name,'volume m3']
        V[seg_name+'_f'] = 0.2 * v_df.loc[seg_name,'volume m3']

    # make a time series of average tracer in the volume defined
    # by the segments of this experiment
    if 'hood_canal' in infile:
        this_seg_list_s = ['H'+ str(item) + '_s' for item in range(3,9)]
        this_seg_list_f = ['H'+ str(item) + '_f' for item in range(3,9)]
        this_seg_list = this_seg_list_s + this_seg_list_f
    
    this_aa = aa.loc[:,this_seg_list]
    this_V = V[this_seg_list]
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

    ax.text(.95,.85, 'e-Folding Time = %d Days' % (td[ind_ef]), fontsize=fs+2,
        transform=ax.transAxes, horizontalalignment='right', color='r', fontweight='bold')

    plt.show()

    fig.savefig(outdir + 'HC_residence_time_plot_' + season + '.png')

