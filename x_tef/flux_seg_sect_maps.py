"""
Plot a map of the sections, segments, and rivers.

"""

# imports
import matplotlib.pyplot as plt
import pickle
import netCDF4 as nc
import pandas as pd
import numpy as np

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun

gridname = 'cas6'
tag = 'v3'
Ldir = Lfun.Lstart(gridname, tag)

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

sys.path.append(os.path.abspath(Ldir['parent'] + 'ptools/pgrid'))
import gfun
import gfun_plotting as gfp
Gr = gfun.gstart(gridname=gridname) # gives access to river info

import tef_fun

import flux_fun
from importlib import reload
reload(flux_fun)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

# select the indir
indir0 = Ldir['LOo'] + 'tef/'
voldir = indir0 + 'volumes_' + Ldir['gridname'] + '/'
# and set the outdir
outdir = indir0 + 'misc_figs_' + Ldir['gridname'] + '/'
Lfun.make_dir(outdir)

# load DataFrame of volumes, created by flux_get_vol.py
v_df = pd.read_pickle(voldir + 'volumes.p')
# index is ['J1', 'J2', 'J3',...
# columns are ['volume m3', 'area m2', 'lon', 'lat']

# load DataFrame of rivers
ri_fn = Gr['ri_dir'] + 'river_info.csv'
riv_df = pd.read_csv(ri_fn, index_col='rname')

def inax(x,y,aa):
    # returns True if x,y, is inside of the box defined by aa
    is_in = x>aa[0] and x<aa[1] and y>aa[2] and y<aa[3]
    return is_in

aa0 = [-125.4, -122, 46.8, 50.4] # Salish Sea
aa1 = [-123.25, -122.1, 47, 48.5] # Puget Sound

plt.close('all')

for do_rivers in [False, True]:
    
    # plotting
    lw=3
    fs=16
    plt.rc('font', size=fs)

    fig = plt.figure(figsize=(16,10))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    
    ax_counter = 0
    for ax in [ax1, ax2]:
    
        if ax_counter == 0:
            aa = aa0
        elif ax_counter == 1:
            aa = aa1

        # sections
        for sn in sect_df.index:
            x0, x1, y0, y1, landward = sect_df.loc[sn,:]
            xx = (x0+x1)/2; yy = (y0+y1)/2
            # color
            for ch in flux_fun.channel_dict.keys():
                if sn in flux_fun.channel_dict[ch]:
                    sn_color = flux_fun.c_dict[ch]
                else:
                    pass
            if inax(xx,yy,aa):
                ax.plot([x0,x1], [y0,y1], '-', color=sn_color, lw=lw)
                if (ax_counter == 0) and not inax(xx,yy,aa1):
                    ax.text(xx,yy, sn, fontsize=.6*fs, color='k', va='center',ha='center',
                    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5), style='italic', weight='bold')
                elif ax_counter == 1:
                    ax.text(xx,yy, sn, fontsize=.6*fs, color='k', va='center',ha='center',
                    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5), style='italic', weight='bold')
    
        # segments
        for Sn in v_df.index:
            vx = v_df.loc[Sn,'lon']
            vy = v_df.loc[Sn,'lat']
            # color
            for ch in flux_fun.short_seg_dict.keys():
                if Sn in flux_fun.short_seg_dict[ch]:
                    Sn_color = flux_fun.c_dict[ch]
                else:
                    pass
            if inax(vx,vy,aa):
                if (ax_counter == 0) and not inax(vx,vy,aa1):
                    # print(ax_counter)
                    # print(Sn)
                    # print(inax(vx,vy,aa1))
                    # print('')
                    if Sn in ['T1', 'T2']:
                        ax.text(vx,vy, Sn,
                            va='center',ha='center', color=Sn_color,
                                fontsize=0.5*fs, fontweight='bold')
                    else:
                        ax.text(vx,vy, Sn,
                            va='center',ha='center', color=Sn_color,
                                fontsize=fs, fontweight='bold')
                        
                elif ax_counter == 1:
                    if Sn in ['T1', 'T2']:
                        ax.text(vx,vy, Sn,
                            va='center',ha='center', color=Sn_color,
                                fontsize=0.5*fs, fontweight='bold')
                    else:
                        ax.text(vx,vy, Sn,
                            va='center',ha='center', color=Sn_color,
                                fontsize=fs, fontweight='bold')
                
        # rivers
        if do_rivers:
            for rn in riv_df.index:
                try:
                    fn_tr = Gr['ri_dir'] + 'tracks/' + rn + '.csv'
                    df_tr = pd.read_csv(fn_tr, index_col='ind')
                    rx = df_tr['lon'].values
                    ry = df_tr['lat'].values
                    if inax(rx[0], ry[0], aa):
                        ax.plot(rx, ry, '-',color='g', lw=lw-1, alpha=.4)
                        ax.plot(rx[-1], ry[-1], '*g', alpha=.4)
                        if (ax_counter == 0) and not inax(rx[-1], ry[-1],aa1):
                            ax.text(rx[-1]+.01, ry[-1]+.01, rn.replace('_',' ').title(),
                                alpha=.4, fontsize=.7*fs, rotation=30)
                        elif ax_counter == 1:
                            ax.text(rx[-1]+.01, ry[-1]+.01, rn.replace('_',' ').title(),
                                alpha=.4, fontsize=.7*fs, rotation=30)
                except FileNotFoundError:
                    pass
        
        pfun.add_coast(ax, color='gray')
        pfun.dar(ax)
        ax.axis(aa)
    
        if ax_counter == 0:
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            ax.set_xticks([-125, -124, -123])
            ax.set_yticks([47, 48, 49, 50])
            pfun.draw_box(ax1, aa1, alpha=.3, linewidth=3)
            ax.set_title('(a) Salish Sea')
            iii = 0
            for ch in flux_fun.channel_list:
                ax.text(.02, .2-iii*.05, ch, color = flux_fun.c_dict[ch],
                transform=ax.transAxes, weight='bold', size=.8*fs)
                iii += 1
        elif ax_counter == 1:
            ax.set_xlabel('Longitude')
            ax.set_xticks([-123, -122.5])
            ax.set_yticks([47, 48])
            ax.set_title('(b) Puget Sound')
    
        ax_counter += 1

    fig.tight_layout()
    plt.show()
    if do_rivers:
        fig.savefig(outdir + 'seg_sect_maps_with_rivers.png')
    else:
        fig.savefig(outdir + 'seg_sect_maps_no_rivers.png')
        
    plt.rcdefaults()

    