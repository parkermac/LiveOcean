"""
Plot the mean of many TEF extractions on all the channels.

"""

# imports
import matplotlib.pyplot as plt
import pickle
import netCDF4 as nc
import pandas as pd
import numpy as np

import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
gridname = 'cas6'; tag = 'v3'
Ldir = Lfun.Lstart(gridname, tag)

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

import tef_fun
import flux_fun
from importlib import reload
reload(flux_fun)

# colors
clist = flux_fun.clist

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

# select input directory
indir0 = Ldir['LOo'] + 'tef/'
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/flux/'

# hacky way of getting the year, assumes "item" is of the form:
# 'cas6_v3_lo8b_2017.01.01_2017.12.31'
year_str = item.split('_')[-1].split('.')[0]
year = int(year_str)

outdir = indir0 + '/tef_all_sections/'
Lfun.make_dir(outdir)

def plotit(ax, sect_df, sect_list, lcol, qsign):
    counter = 0
    for sn in sect_list:
        # some information about direction
        x0, x1, y0, y1, landward = sect_df.loc[sn,:]
        ax.plot([x0,x1], [y0,y1], '-', color=lcol, linewidth=3)
        xx = (x0+x1)/2
        yy = (y0+y1)/2
        # add a mark showing which direction the deep flow is going,
        # noting that the transports are now all positive east or north.
        clat = np.cos(np.pi*yy/180)
        if (x0==x1) and (y0!=y1):
            sdir = 'NS'
            dd = qsign[counter] * 0.05 / clat
            ww = dd/4
            ax.fill([xx, xx+dd, xx], [yy-ww, yy, yy+ww], color=lcol)
        elif (x0!=x1) and (y0==y1):
            sdir = 'EW'
            dd = qsign[counter] * 0.05
            ww = dd/(4*clat)
            ax.fill([xx-ww, xx, xx+ww], [yy, yy+dd, yy], color=lcol)
        counter += 1

# PLOTTING
plt.close('all')
channel_list = flux_fun.channel_list
channel_dict = flux_fun.long_channel_dict
lcol_dict = dict(zip(flux_fun.channel_list, clist))
channel_list.reverse() # makes overlaying colors look better
lw=2
fs=16
distmax = 420

for season in flux_fun.season_list:
    
    # load the "season" DataFrame made by flux_make_two_layer.py
    # which has columns=['q_s', 'q_f', 'f_s', 'f_f', 's_s', 's_f', 'lon', 'lat']
    # and the index is all the section names
    fn = indir + 'two_layer_' + season + '.p'
    df = pickle.load(open(fn,'rb'))

    fig = plt.figure(figsize=(15,10))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(122) # map


    # create all the distance vectors and save in a dict
    dist_dict = {}
    for ch_str in channel_list:
        sect_list = channel_dict[ch_str]
        x = df.loc[sect_list,'lon'].to_numpy(dtype='float')
        y = df.loc[sect_list,'lat'].to_numpy(dtype='float')
        dist_dict[ch_str] = flux_fun.make_dist(x,y)
    
    # adjust the distance vectors to join at the correct locations
    ind_ai = channel_dict['Juan de Fuca to Strait of Georgia'].index('jdf4')
    dist0_ai = dist_dict['Juan de Fuca to Strait of Georgia'][ind_ai]
    dist_dict['Admiralty Inlet to South Sound'] += dist0_ai
    #
    ind_hc = channel_dict['Admiralty Inlet to South Sound'].index('ai3')
    dist0_hc = dist_dict['Admiralty Inlet to South Sound'][ind_hc]
    dist_dict['Hood Canal'] += dist0_hc
    #
    ind_wb = channel_dict['Admiralty Inlet to South Sound'].index('ai4')
    dist0_wb = dist_dict['Admiralty Inlet to South Sound'][ind_wb]
    dist_dict['Whidbey Basin'] += dist0_wb
 
    do_plot_extras = True
    for ch_str in channel_list:
        print('== ' + ch_str + ' ==')
        sect_list = channel_dict[ch_str]
    
        q_s = df.loc[sect_list,'q_s'].to_numpy(dtype='float')/1e3
        q_f = df.loc[sect_list,'q_f'].to_numpy(dtype='float')/1e3
        s_s = df.loc[sect_list,'s_s'].to_numpy(dtype='float')
        s_f = df.loc[sect_list,'s_f'].to_numpy(dtype='float')
    
        dist = dist_dict[ch_str]
    
        lcol = lcol_dict[ch_str]
        ax1.plot(dist,np.abs(q_s),'-', color=lcol,linewidth=lw, label=ch_str)
        ax1.plot(dist,np.abs(q_f),'-', color=lcol,linewidth=lw, label=ch_str)
        ax1.set_xlim(-5,distmax)
        ax1.set_ylim(0, 200)
        ax1.grid(True)
        ax1.set_ylabel('Qin and Qout (1000 m3/s)', fontsize=fs)
        year_str = indir.split('/')[-3].split('_')[-1].split('.')[0] # brittle!
        ax1.set_title(season.title() + ' ' + year_str, fontsize=fs)

        counter = 0
        for sn in sect_list:
            sn = sn.upper()
            ax1.text(dist[counter], np.abs(q_s[counter]), sn, rotation=45, fontsize=8)
            counter += 1
        
        ax2.fill_between(dist, s_s, y2=s_f, color=lcol, alpha=.5)
        ax2.set_xlim(-5,distmax)
        ax2.set_ylim(21.5,34)
        ax2.grid(True)
        ax2.set_xlabel('Distance from Mouth (km)', fontsize=fs)
        ax2.set_ylabel('Salinity', fontsize=fs)
    
        if do_plot_extras:
            aa = [-125.5, -122, 46.7, 50.4]
            pfun.add_coast(ax3)
            pfun.dar(ax3)
            ax3.axis(aa)
            ax3.set_title('Section Locations, and direction of deep inflow')
        
            for sn in sect_df.index:
                x0, x1, y0, y1, landward = sect_df.loc[sn,:]
                xx = (x0+x1)/2
                yy = (y0+y1)/2
                ax3.text(xx,yy, sn, rotation=45, fontsize=8)
            do_plot_extras = False
        
        # get the sign of q_s and plot the section locations with "inflow" direction
        # (defined as the direction of transport of the saltier water of the pair)
        qsign = np.sign(q_s)
        plotit(ax3, sect_df, sect_list, lcol, qsign)

    #plt.show()
    fig.savefig(outdir + 'all_sections_' + year_str + '_' + season + '.png')
    