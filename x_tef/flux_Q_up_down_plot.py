"""
Code to plot the spatial distribution of the vertical velocities that
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

from importlib import reload
reload(flux_fun)

# Input directory
indir0 = Ldir['LOo'] + 'tef/'
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/flux/'

outdir = indir0 + item + '/misc_figs/'

# load a Series of the volumes of each segment, created by flux_get_vol.py
v_df = pd.read_pickle(indir + 'volumes.p')
# index is ['J1', 'J2', 'J3',...
# columns are ['volume m3', 'area m2', 'lon', 'lat']

clist = flux_fun.clist
channel_list = flux_fun.channel_list
channel_dict = flux_fun.channel_dict
lcol_dict = dict(zip(flux_fun.channel_list, clist))

seg_dict = flux_fun.short_seg_dict

plt.close('all')

for season in flux_fun.dtr.keys():

    df_dict = {}
    #for ch_str in ['Whidbey Basin']:
    #for ch_str in ['Juan de Fuca to Strait of Georgia']:
    for ch_str in channel_list:

            # this is the big DataFrame created by flux_get_A.py
            q_seg_df = pd.read_pickle(indir + 'q_df_' + season + '.p')
            # index is ['J1_s', 'J1_f', 'J2_s',... = (*)
            # columns are ['ocean_s', 'ocean_f', 'river_s', 'river_f', 'J1_s', 'J1_f', 'J2_s',...
    
            # load the "season" DataFrame made by flux_make_two_layer.py
            # which has columns=['q_s', 'q_f', 'f_s', 'f_f', 's_s', 's_f', 'lon', 'lat']
            # and the index is all the section names
            fn = indir + 'two_layer_' + season + '.p'
            q_sect_df = pickle.load(open(fn,'rb'))
    
            # working it out channel by channel
            # In this case we start on a segment and end on a section
            # and the two lists are the same length
            sect_list = channel_dict[ch_str].copy()
            seg_list = seg_dict[ch_str].copy()
            NSC = len(sect_list)
            NSG = len(seg_list)
    
            x = q_sect_df.loc[sect_list,'lon'].to_numpy(dtype='float')
            y = q_sect_df.loc[sect_list,'lat'].to_numpy(dtype='float')
            dist = flux_fun.make_dist(x,y)
    
            # make vectors of net vertical transport on the sections of this channel
            # NSS = len(seg_list)
            q_s = np.nan * np.ones(NSC)
            q_f = np.nan * np.ones(NSC)
            s_s = np.nan * np.ones(NSC)
            s_f = np.nan * np.ones(NSC)

            qnet_up = np.nan * np.ones(NSC)
            qnet_down = np.nan * np.ones(NSC)
            qnet_r = np.nan * np.ones(NSC)
            qnet = np.nan * np.ones(NSC)
            
            q_up = np.nan * np.ones(NSG)
            q_down = np.nan * np.ones(NSG)
            q_r = np.nan * np.ones(NSG)

            counter = 0
            for sect in sect_list:
                q_s[counter] = q_sect_df.loc[sect,'q_s']
                q_f[counter] = q_sect_df.loc[sect,'q_f']
                s_s[counter] = q_sect_df.loc[sect,'s_s']
                s_f[counter] = q_sect_df.loc[sect,'s_f']
                counter += 1
                
            counter = 0
            for seg in seg_list:
                q_up[counter] = q_seg_df.loc[seg+'_f',seg+'_s']
                q_down[counter] = q_seg_df.loc[seg+'_s',seg+'_f']
                q_r[counter] = q_seg_df.loc[seg+'_f','river_f']
                counter += 1
            
            q_s = np.abs(q_s)/1e3
            q_f = np.abs(q_f)/1e3
            q_up = q_up/1e3
            q_down = q_down/1e3
            q_r = q_r/1e3
            
            if NSC > NSG:
                print(50*'=')
                print('Channel with two open ends')
                print(ch_str)
                
                # "_a" = sum from "landward direction"
                # but still the first point is on the ocean end
                qnet_up_a = np.cumsum(q_up[::-1])[::-1]
                qnet_down_a = np.cumsum(q_down[::-1])[::-1]
                qnet_r_a = np.cumsum(q_r[::-1])[::-1]
                qnet_a = qnet_up_a - qnet_down_a
                # "_b" = sum from seaward direction
                qnet_up_b = np.cumsum(q_up)
                qnet_down_b = np.cumsum(q_down)
                qnet_r_b = np.cumsum(q_r)
                qnet_b = qnet_up_b - qnet_down_b

                dQ_a = qnet_a[0] - q_s[0]
                dQ_b = qnet_b[-1] - q_s[-1]
                
                qnet_aa = np.append(qnet_a, 0)
                qnet_bb = np.append(0, qnet_b)
                
                qnet_aa -= dQ_a
                qnet_bb -= dQ_b
                
                # now let's find a way to interpolate to get the amount to split
                i0, i1, fr = zfun.get_interpolant(np.array([0]), qnet_bb)
                i0 = i0[0]
                i1 = i1[0]
                print('i0=%d, i1=%d, fr=%0.2f' % (i0, i1, fr))
                
                # qnet[0:i0+1] = qnet_aa[0:i0+1]
                # qnet[i1:] = qnet_bb[i1:]
                qnet[0:i0+1] = qnet_aa[0:i0+1] - (qnet_aa[i0] + fr*(qnet_aa[i1]-qnet_aa[i0]))
                qnet[i1:] = qnet_bb[i1:] - (qnet_bb[i0] + fr*(qnet_bb[i1]-qnet_bb[i0]))
                
                # Now apply the same split to the other vectors
                qnet_up_aa = np.append(qnet_up_a, 0)
                qnet_up_bb = np.append(0, qnet_up_b)
                qnet_down_aa = np.append(qnet_down_a, 0)
                qnet_down_bb = np.append(0, qnet_down_b)
                qnet_r_aa = np.append(qnet_r_a, 0)
                qnet_r_bb = np.append(0, qnet_r_b)
                
                qnet_up[0:i0+1] = qnet_up_aa[0:i0+1] - (qnet_up_aa[i0] + fr*(qnet_up_aa[i1]-qnet_up_aa[i0]))
                qnet_up[i1:] = qnet_up_bb[i1:] - (qnet_up_bb[i0] + fr*(qnet_up_bb[i1]-qnet_up_bb[i0]))

                qnet_down[0:i0+1] = qnet_down_aa[0:i0+1] - (qnet_down_aa[i0] + fr*(qnet_down_aa[i1]-qnet_down_aa[i0]))
                qnet_down[i1:] = qnet_down_bb[i1:] - (qnet_down_bb[i0] + fr*(qnet_down_bb[i1]-qnet_down_bb[i0]))

                qnet_r[0:i0+1] = qnet_r_aa[0:i0+1] - (qnet_r_aa[i0] + fr*(qnet_r_aa[i1]-qnet_r_aa[i0]))
                qnet_r[i1:] = qnet_r_bb[i1:] - (qnet_r_bb[i0] + fr*(qnet_r_bb[i1]-qnet_r_bb[i0]))
                
            else:
                
                qnet_up = np.cumsum(q_up[::-1])[::-1]
                qnet_down = np.cumsum(q_down[::-1])[::-1]
                qnet_r = np.cumsum(q_r[::-1])[::-1]
                qnet = qnet_up - qnet_down
        
            df = pd.DataFrame(index=sect_list)
            df.loc[:,'dist'] = dist
            df.loc[:,'q_s'] = q_s
            df.loc[:,'q_f'] = q_f
            df.loc[:,'s_s'] = s_s
            df.loc[:,'s_f'] = s_f
            df.loc[:,'qnet_up'] = qnet_up
            df.loc[:,'qnet_down'] = qnet_down
            df.loc[:,'qnet'] = qnet
            df.loc[:,'qnet_r'] = qnet_r
        
            df_dict[ch_str] = df
        
    # connections
    for qn in ['qnet_up', 'qnet_down', 'qnet', 'qnet_r']:
        df_dict['Admiralty Inlet to South Sound'].loc[['ai1','ai2','ai3','ai4'],qn] += df_dict['Whidbey Basin'].loc['wb1',qn]
        df_dict['Admiralty Inlet to South Sound'].loc[['ai1','ai2','ai3'],qn] += df_dict['Hood Canal'].loc['hc1',qn]
        df_dict['Juan de Fuca to Strait of Georgia'].loc[['jdf1','jdf2','jdf3','jdf4'],qn] += df_dict['Admiralty Inlet to South Sound'].loc['ai1',qn]


    qlim_dict = {'Juan de Fuca to Strait of Georgia':400,
     'Admiralty Inlet to South Sound':100,
     'Hood Canal':10,
     'Whidbey Basin':10}

    fig = plt.figure(figsize=(11,8))
    fig.suptitle(season.title())
    ii = 1
    for ch_str in channel_dict:

        df = df_dict[ch_str]

        ax = fig.add_subplot(3,4,ii)
        
        if ii==1:
            legend=True
        else:
            legend=False
        df.plot(x='dist', y=['q_s', 'q_f', 'qnet_up', 'qnet_down'],#, 'qnet'],
            grid=True, ax=ax, ylim=(0,qlim_dict[ch_str]), legend=legend)
        ax.set_title(ch_str, fontsize=8)
        for sn in df.index:
            ax.text(df.loc[sn,'dist'],df.loc[sn,'q_s'],sn, fontsize=6, alpha=.5,
            horizontalalignment='center')
        ax.set_xlabel('')
        ax.set_xticklabels([])

        ax = fig.add_subplot(3,4,ii+4)
        df.plot(x='dist', y=['s_s', 's_f'], grid=True, ax=ax, ylim=(22,34), legend=legend)
        ax.set_xlabel('')
        ax.set_xticklabels([])

        ax = fig.add_subplot(3,4,ii+8)
        df.plot(x='dist', y=['qnet_r'], grid=True, ax=ax, ylim=(0,8), legend=legend)

        ii += 1

    plt.show()

    fig.savefig(outdir + 'Q_up_down_' + season + '.png')

    
