"""
Plot the results of the flux engine.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pickle
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

# Output directory
outdir = indir0 + 'validation_plots/'
Lfun.make_dir(outdir)

plt.close('all')

fs = 16
plt.rc('font', size=fs)

for year in [2017, 2018, 2019]:
    year_str = str(year)
    item = 'cas6_v3_lo8b_'+year_str+'.01.01_'+year_str+'.12.31'
    indir = indir0 + item + '/flux/'
    indir2 = indir0 + 'flux_engine/'
    
    season = 'full'
    # this is only valid for the full year - seasons are not in equilibrium
    
    # load the DataFrame of TEF transports and salinities, created by flux_make_two_layer.py
    df_2 = pd.read_pickle(indir + 'two_layer_' + season + '.p')
    # index is ['jdf1', 'jdf2', 'jdf3',...
    # columns are ['q_s', 'q_f', 'f_s', 'f_f', 's_s', 's_f', 'lon', 'lat']

    # load the segment definitions, which is a dictionary defining all the
    # sections and rivers surrounding each segment
    segs = flux_fun.segs
    # this is a dict of dicts, with entries like:
    #   {'J1':
    #       {'S': [], 'N': [], 'W': ['jdf1'], 'E': ['jdf2'], 'R': ['sanjuan', 'hoko']},...

    # load the DataFrame of results of flux_engine.py
    cc = pd.read_pickle(indir2 + 'cas6_v3_lo8b/S_OceanSalt_' + str(year) +'_'+ season + '_AGE.p')
    # index is ['J1_s', 'J1_f',... = (*)
    # columns are ['c', 'v', 'netq']

    fig = plt.figure(figsize=(13,8))

    ax_counter = 1
    for ch in flux_fun.channel_dict.keys():

        if ax_counter == 1:
            ax = fig.add_subplot(2,1,1)
            ax.set_xlim(-10,410)
        elif ax_counter == 2:
            ax = fig.add_subplot(2,2,3)
            ax.set_xlim(-10,170)
        elif ax_counter == 3:
            ax = fig.add_subplot(2,4,7)
            ax.set_xlim(-10,100)
        elif ax_counter == 4:
            ax = fig.add_subplot(2,4,8)
            ax.set_xlim(-10,100)
    
        sect_list = flux_fun.channel_dict[ch]
        seg_list = flux_fun.seg_dict[ch]

        # get target salinities (the actual TEF values, at segment boundaries)
        Ns = len(sect_list)
        sa_s = np.zeros(Ns)
        sa_f = np.zeros(Ns)
        sa_x = np.zeros(Ns)
        sa_y = np.zeros(Ns)
        counter = 0
        for sect in sect_list:
            q_s, q_f, f_s, f_f, s_s, s_f, lon, lat = df_2.loc[sect,:]
            sa_s[counter] = s_s
            sa_f[counter] = s_f
            sa_x[counter] = lon
            sa_y[counter] = lat
            counter += 1
        # and associated distance vector
        sa_dist = flux_fun.make_dist(sa_x, sa_y)

        vs = [s + '_s' for s in seg_list]
        vf = [s + '_f' for s in seg_list]

        # alternate definition of dist that goes with segments
        if ch in ['Admiralty Inlet to South Sound', 'Hood Canal']:
            dist = sa_dist.copy()
            ddd = np.diff(sa_dist)/2
            dist[:-1] += ddd
            dist[-1] += ddd[-1]
            dist = np.append(-5,dist)
        elif ch == 'Whidbey Basin':
            dist = sa_dist[:-1].copy()
            ddd = np.diff(sa_dist)/2
            dist += ddd
            dist = np.append(-5,dist)
        else:
            dist = sa_dist[:-1].copy()
            ddd = np.diff(sa_dist)/2
            dist += ddd
        
        # Plots of layer salinities

        # values from the flux_engine
        ax.fill_between(dist,cc.loc[vs,'c'],cc.loc[vf,'c'],color='b',alpha=.5)
        
        # TEF target values
        ax.fill_between(sa_dist,sa_s,sa_f,color='r',alpha=.5)

        if False:
            for ii in range(len(dist)):
                ax.text(dist[ii], cc.loc[vs[ii],'c'], seg_list[ii], color='b')
            for ii in range(len(sa_dist)):
                ax.text(sa_dist[ii], sa_s[ii], sect_list[ii], color='r')
            
        ax.set_ylim(22,34)
        
        if ax_counter in [2,3,4]:
            ax.set_xlabel('Distance [km]')
            
        if ax_counter in [3,4]:
            ax.set_yticklabels([])
            
        if ax_counter in [1,2]:
            ax.set_ylabel('Salinity')
            
        ax.set_yticks([24,28,32])
            
        if ax_counter == 1:
            ax.text(.9,.55,'TEF Sections',color='r',fontweight='bold',
                transform=ax.transAxes, horizontalalignment='right',alpha=.5)
            ax.text(.9,.75,'Box Model Segments',color='b',fontweight='bold',
                transform=ax.transAxes, horizontalalignment='right',alpha=.5)
            ax.text(.9,.1,year_str,color='k',fontweight='bold',style='italic',
                transform=ax.transAxes, horizontalalignment='right')
                
        abc = 'abcd'
        ax.text(.05,.05,'(%s) %s' % (abc[ax_counter-1],ch),
            transform=ax.transAxes, color=flux_fun.c_dict[ch],weight='bold')

        ax.grid(True)

        ax_counter += 1
        
    fig.tight_layout()
    fig.savefig(outdir + 'validation_plot_' + year_str + '_' + season + '.png')

plt.show()
plt.rcdefaults()


