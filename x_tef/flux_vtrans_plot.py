"""
Code to plot the spatial distribution of the vertical transports that
result from flux_get_A.py.

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
outdir = indir0 + 'vertical_transport_plots/'
Lfun.make_dir(outdir)

voldir = indir0 + 'volumes_' + Ldir['gridname'] + '/'
# load a Series of the volumes of each segment, created by flux_get_vol.py
v_df = pd.read_pickle(voldir + 'volumes.p')
# index is ['J1', 'J2', 'J3',...
# columns are ['volume m3', 'area m2', 'lon', 'lat']

plt.close('all')

fs = 16
lw = 3
plt.rc('font', size=fs)

for year in [2017, 2018, 2019]:
    year_str = str(year)
    item = 'cas6_v3_lo8b_'+year_str+'.01.01_'+year_str+'.12.31'
    indir = indir0 + item + '/flux/'
    indir2 = indir0 + 'flux_engine/'
    
    season = 'full'
    # this is only valid for the full year - seasons are not in equilibrium
    
    # this is the big DataFrame created by flux_get_A.py
    q_df = pd.read_pickle(indir + 'q_df_' + season + '.p')
    # index is ['J1_s', 'J1_f', 'J2_s',... = (*)
    # columns are ['ocean_s', 'ocean_f', 'river_s', 'river_f', 'J1_s', 'J1_f', 'J2_s',...
    
    fig = plt.figure(figsize=(14,8))

    ax_counter = 1
    for ch in flux_fun.seg_dict.keys():
    
        if ax_counter == 1:
            ax = fig.add_subplot(2,1,1)
            ax.set_xlim(-5,410)
        elif ax_counter == 2:
            ax = fig.add_subplot(2,2,3)
            ax.set_xlim(-5,200)
        elif ax_counter == 3:
            ax = fig.add_subplot(2,4,7)
            ax.set_xlim(-5,105)
        elif ax_counter == 4:
            ax = fig.add_subplot(2,4,8)
            ax.set_xlim(-5,80)

        sect_list = flux_fun.channel_dict[ch]
        seg_list = flux_fun.short_seg_dict[ch]
        #print(seg_list)
                
        dist = flux_fun.make_dist(v_df.loc[seg_list,'lon'],v_df.loc[seg_list,'lat'])
        
        # make vectors of vertical transport on the segments of this channel
        q_up = np.nan * np.ones(len(seg_list))
        q_down = np.nan * np.ones(len(seg_list))
        w_up = np.nan * np.ones(len(seg_list))
        w_down = np.nan * np.ones(len(seg_list))
        seg_counter = 0
        for seg in seg_list:
            q_up[seg_counter] = q_df.loc[seg+'_f',seg+'_s']
            q_down[seg_counter] = q_df.loc[seg+'_s',seg+'_f']
            w_up[seg_counter] = 86400 *q_df.loc[seg+'_f',seg+'_s'] / v_df.loc[seg,'area m2']
            w_down[seg_counter] = 86400 * q_df.loc[seg+'_s',seg+'_f'] / v_df.loc[seg,'area m2']
            seg_counter += 1 
    
        # plotting
        upcol = 'lightsalmon'
        dncol = 'mediumslateblue'
        yld = {1:125, 2:25, 3:5, 4:5}
        #yld = {1:5, 2:50, 3:5, 4:5} # for vertical velocity
        ax.plot(dist, q_up/1e3,'-o', color=upcol, linewidth=lw)
        ax.plot(dist, -q_down/1e3,'-o', color=dncol, linewidth=lw)
        ax.text(.02,.03,'Net: %0.1f up, %0.1f down' % ((q_up/1e3).sum(),(-q_down/1e3).sum()),
            transform = ax.transAxes, size=.7*fs, weight='bold')
        # ax.plot(dist, w_up,'-o', color=upcol, linewidth=lw)
        # ax.plot(dist, -w_down,'-o', color=dncol, linewidth=lw)
        for ii in range(len(dist)):
            ax.text(dist[ii], -.6*yld[ax_counter], seg_list[ii],
            ha='center',size=.7*fs, alpha=.5)
            
        if ax_counter in [2,3,4]:
            ax.set_xlabel('Distance [km]')
            
        if ax_counter in [3,4]:
            ax.set_yticks([-10,-5,0,5,10])
            
        if ax_counter ==1:
            ax.text(.05,.6,'Upward',color=upcol,fontweight='bold',
                transform=ax.transAxes, va='center')
            ax.text(.05,.4,'Downward',color=dncol,fontweight='bold',
                transform=ax.transAxes, va='center')
            ax.text(.9,.1,year_str,color='k',fontweight='bold',style='italic',
                transform=ax.transAxes, ha='right')
        
        if ax_counter in [1,2]:
            ax.set_ylabel('Vert. Trans. $[10^{3}m^{3}s^{-1}]$')
        
        abc = 'abcd'
        ax.text(.02,.9,'(%s) %s' % (abc[ax_counter-1],ch),
            transform=ax.transAxes, color=flux_fun.c_dict[ch], weight='bold')
        
        ax.axhline(color='k')
        
        ax.set_ylim(-yld[ax_counter],yld[ax_counter])
    
        ax_counter += 1

    fig.tight_layout()

    fig.savefig(outdir + 'vtrans_' + year_str + '.png')

plt.show()
plt.rcdefaults()


