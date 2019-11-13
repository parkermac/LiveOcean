"""
Plot the results of the flux age engine.

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

# select the indir
indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

# load the DataFrame of results of flux_engine.py
infile = Lfun.choose_item(indir, tag='cc_', itext='Choose flux engine output file:')
cc = pd.read_pickle(indir + infile)
cc.loc[:,'age'] = 365*cc.loc[:,'ca']/cc.loc[:,'c']

# load a Series of the volumes of each segment, created by flux_get_vol.py
v_df = pd.read_pickle(indir + 'volumes.p')


#plt.close('all')
fig = plt.figure(figsize=(13,8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

fs = 16
lw = 3

Dist_dict = {}
for ch in flux_fun.channel_dict.keys():
    seg_list = flux_fun.seg_dict[ch]
    dist = flux_fun.make_dist(v_df.loc[seg_list,'lon'],v_df.loc[seg_list,'lat'])
    Dist = pd.Series(dist, index=seg_list)
    Dist_dict[ch] = Dist
# hand edits
Dist_dict['Admiralty Inlet to South Sound'] += Dist_dict['Juan de Fuca to Strait of Georgia']['J4']
Dist_dict['Hood Canal'] += Dist_dict['Admiralty Inlet to South Sound']['A3']
Dist_dict['Whidbey Basin'] += Dist_dict['Admiralty Inlet to South Sound']['M1']

ch_counter = 0
for ch in flux_fun.channel_dict.keys():
    dist = Dist_dict[ch]
    seg_list = flux_fun.seg_dict[ch]

    vs = [s + '_s' for s in seg_list]
    vf = [s + '_f' for s in seg_list]

    color = flux_fun.clist[ch_counter]
    # ax1.plot(dist, cc.loc[vs,'c'].values,'-',color=color, linewidth=lw)
    # ax1.plot(dist, cc.loc[vf,'c'].values,'-',color=color, linewidth=lw, alpha=.5)
    ax1.plot(dist, np.log10(cc.loc[vs,'c'].values),'-',color=color, linewidth=lw)
    ax1.plot(dist, np.log10(cc.loc[vf,'c'].values),'--',color=color, linewidth=lw)
    # for ii in range(len(dist)):
    #     ax1.text(dist[ii],np.log10(cc.loc[vs[ii],'c']), seg_list[ii],color=color)
        
    ax2.plot(dist, cc.loc[vs,'age'].values,'-',color=color, linewidth=lw)
    ax2.plot(dist, cc.loc[vf,'age'].values,'--',color=color, linewidth=lw)
    for ii in range(len(dist)):
        ax2.text(dist[ii], cc.loc[vs[ii],'age'], seg_list[ii], color=color,
        horizontalalignment='center', verticalalignment='bottom', fontsize=fs-3)
        
    ax2.text(.05, .75-ch_counter*.11, ch, color=color,
        transform=ax2.transAxes, fontsize=fs,
        bbox=dict(facecolor='w',edgecolor='w', alpha=0.8))
    
    ch_counter += 1
    
#ax1.set_ylim(bottom=0)
ax2.set_ylim(0, 700)

ax1.set_xlim(-10,410)
ax2.set_xlim(-10,410)

ax1.grid(True)
ax2.grid(True)

ax1.set_ylabel('log10(Concentration)', fontsize=fs)
ax2.set_ylabel('Age (days)', fontsize=fs)
ax2.set_xlabel('Distance from Mouth of JdF (km)', fontsize=fs)

ax1.tick_params(labelsize=fs-2) 
ax2.tick_params(labelsize=fs-2)

ttext = infile.replace('cc_','').replace('.p','').replace('_',' ').title()
ax1.set_title('Source: ' + ttext, fontsize=fs)

plt.show()


