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

for season in ['full']:#flux_fun.dtr.keys():

    df_dict = {}
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
            sect_list.reverse()
            seg_list.reverse()
    
            if ch_str in ['Juan de Fuca to Strait of Georgia', 'Whidbey Basin']:
                sect_list = sect_list[1:]
        
            ss_dict = dict(zip(sect_list,seg_list))
    
            x = q_sect_df.loc[sect_list,'lon'].to_numpy(dtype='float')
            y = q_sect_df.loc[sect_list,'lat'].to_numpy(dtype='float')
            dist = flux_fun.make_dist(x,y)
    
            # make vectors of net vertical transport on the segments of this channel
            NSS = len(seg_list)
            q_up = np.nan * np.ones(NSS)
            q_down = np.nan * np.ones(NSS)
            q_s = np.nan * np.ones(NSS)
            q_f = np.nan * np.ones(NSS)
            s_s = np.nan * np.ones(NSS)
            s_f = np.nan * np.ones(NSS)
            q_r = np.nan * np.ones(NSS)

            counter = 0
            for sect in sect_list:
                seg = ss_dict[sect]
                q_s[counter] = q_sect_df.loc[sect,'q_s']
                q_f[counter] = q_sect_df.loc[sect,'q_f']
                s_s[counter] = q_sect_df.loc[sect,'s_s']
                s_f[counter] = q_sect_df.loc[sect,'s_f']
                q_up[counter] = q_seg_df.loc[seg+'_f',seg+'_s']
                q_down[counter] = q_seg_df.loc[seg+'_s',seg+'_f']
                q_r[counter] = q_seg_df.loc[seg+'_f','river_f']
                counter += 1
        
            qnet_up = np.cumsum(q_up)
            qnet_down = np.cumsum(q_down)
            qnet_r = np.cumsum(q_r)
        
            df = pd.DataFrame(index=seg_list)
            df.loc[:,'sect_list'] = sect_list
            df.loc[:,'dist'] = dist
            df.loc[:,'q_s'] = np.abs(q_s)/1e3
            df.loc[:,'q_f'] = np.abs(q_f)/1e3
            df.loc[:,'s_s'] = s_s
            df.loc[:,'s_f'] = s_f
            df.loc[:,'qnet_up'] = qnet_up/1e3
            df.loc[:,'qnet_down'] = qnet_down/1e3
            df.loc[:,'qnet'] = (qnet_up-qnet_down)/1e3
            df.loc[:,'qnet_r'] = qnet_r/1e3
        
            df_dict[ch_str] = df
        
    # connections
    for qn in ['qnet_up', 'qnet_down', 'qnet', 'qnet_r']:
        df_dict['Admiralty Inlet to South Sound'].loc[['M1','A3','A2','A1'],qn] += df_dict['Whidbey Basin'].loc['W1',qn]
        df_dict['Admiralty Inlet to South Sound'].loc[['A3','A2','A1'],qn] += df_dict['Hood Canal'].loc['H1',qn]
        df_dict['Juan de Fuca to Strait of Georgia'].loc[['J4','J3','J2','J1'],qn] += df_dict['Admiralty Inlet to South Sound'].loc['A1',qn]
        
    # adjustments to account for switched flow directions
    df_dict['Juan de Fuca to Strait of Georgia'].loc[:,'qnet_up'] += df_dict['Juan de Fuca to Strait of Georgia'].loc['G5','qnet']
    df_dict['Juan de Fuca to Strait of Georgia'].loc[:,'qnet_down'] += df_dict['Juan de Fuca to Strait of Georgia'].loc['G5','qnet']
    df_dict['Juan de Fuca to Strait of Georgia'].loc[:,'qnet_up'] = df_dict['Juan de Fuca to Strait of Georgia'].loc[:,'qnet_up'].abs()
    df_dict['Juan de Fuca to Strait of Georgia'].loc[:,'qnet_down'] = df_dict['Juan de Fuca to Strait of Georgia'].loc[:,'qnet_down'].abs()
    df_dict['Juan de Fuca to Strait of Georgia'].loc[:,'qnet'] = (
        df_dict['Juan de Fuca to Strait of Georgia'].loc[:,'qnet_up'] - df_dict['Juan de Fuca to Strait of Georgia'].loc[:,'qnet_down'])

    qlim_dict = {'Juan de Fuca to Strait of Georgia':400,
     'Admiralty Inlet to South Sound':100,
     'Hood Canal':10,
     'Whidbey Basin':10}
 
    fig = plt.figure(figsize=(14,14))
    fig.suptitle(season.title())
    ii = 1
    for ch_str in channel_dict:
    
        df = df_dict[ch_str]
    
        ax = fig.add_subplot(3,4,ii)
        df.plot(x='dist', y=['q_s', 'q_f', 'qnet_up', 'qnet_down', 'qnet'], grid=True, title=ch_str, ax=ax, ylim=(0,qlim_dict[ch_str]))

        ax = fig.add_subplot(3,4,ii+4)
        df.plot(x='dist', y=['s_s', 's_f'], grid=True, ax=ax, ylim=(22,34))
    
        ax = fig.add_subplot(3,4,ii+8)
        df.plot(x='dist', y=['qnet_r'], grid=True, ax=ax, ylim=(0,8))
    
        ii += 1
    
    plt.show()
    
    fig.savefig(outdir + 'Q_up_down_' + season + '.png')

    
