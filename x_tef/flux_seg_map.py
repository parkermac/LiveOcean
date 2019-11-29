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
from importlib import reload
reload(tef_fun)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

# select the indir
indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

outdir = indir0 + 'misc_figs/'
Lfun.make_dir(outdir)

# load DataFrame of volumes, created by flux_get_vol.py
v_df = pd.read_pickle(indir + 'volumes.p')
# index is ['J1', 'J2', 'J3',...
# columns are ['volume m3', 'area m2', 'lon', 'lat']

# load DataFrame of rivers
ri_fn = Gr['ri_dir'] + 'river_info.csv'
riv_df = pd.read_csv(ri_fn, index_col='rname')

# plotting
plt.close('all')
lw=3
fs=14
fig1 = plt.figure(figsize=(9,12))
fig2 = plt.figure(figsize=(9,12))
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)

def inax(x,y,aa, aa_x):
    # returns true if x,y, is inside of the box defined by aa
    is_in = x>aa[0] and x<aa[1] and y>aa[2] and y<aa[3]
    if len(aa_x) == 4:
        is_in_x = x>aa_x[0] and x<aa_x[1] and y>aa_x[2] and y<aa_x[3]
    else:
        is_in_x = False
        
    return is_in and not is_in_x

aa0 = [-125.5, -122, 46.7, 50.4]
aa1 = [-123.3, -122, 46.9, 48.5]

ax_counter = 0
for ax in [ax1, ax2]:
    
    
    if ax_counter == 0:
        aa = aa0
        aa_x = []
    elif ax_counter == 1:
        aa = aa1
        aa_x = []

    # sections
    for sn in sect_df.index:
        x0, x1, y0, y1, landward = sect_df.loc[sn,:]
        xx = (x0+x1)/2; yy = (y0+y1)/2
        if inax(xx,yy,aa, aa_x):
            ax.plot([x0,x1], [y0,y1], '-', color='purple', linewidth=lw)
            ax.text(xx,yy, sn, rotation=45, fontsize=fs-4, color='purple')
    
    # segments
    for Sn in v_df.index:
        vx = v_df.loc[Sn,'lon']
        vy = v_df.loc[Sn,'lat']
        if inax(vx,vy,aa, aa_x):
            ax.text(vx,vy, Sn,
                verticalalignment='center',horizontalalignment='center',
                    fontsize=fs, fontweight='bold')
                
    # rivers
    for rn in riv_df.index:
        try:
            fn_tr = Gr['ri_dir'] + 'tracks/' + rn + '.csv'
            df_tr = pd.read_csv(fn_tr, index_col='ind')
            rx = df_tr['lon'].values
            ry = df_tr['lat'].values
            if inax(rx[0], ry[0], aa, aa_x):
                ax.plot(rx, ry, '-',color='g', linewidth=lw-1, alpha=.4)
                ax.plot(rx[-1], ry[-1], '*g', alpha=.4)
                ax.text(rx[-1]+.01, ry[-1]+.01, rn.replace('_',' ').title(),
                    alpha=.4, fontsize=fs-4, rotation=30)
        except FileNotFoundError:
            pass
        
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    
    ax_counter += 1


plt.show()

fig1.savefig(outdir + 'seg_map_full.png')
fig2.savefig(outdir + 'seg_map_focus.png')
    