"""
Plot the mean of many TEF extractions on a thalweg section.
"""

# imports
import matplotlib.pyplot as plt
import pickle
import netCDF4 as nc
import pandas as pd

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart('cas4', 'v2')

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

pth = os.path.abspath(Ldir['parent'] + 'ptools/pgrid')
if pth not in sys.path:
    sys.path.append(pth)
import gfun
import gfun_plotting as gfp
Gr = gfun.gstart('cas4')

import tef_fun
from importlib import reload
reload(tef_fun)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

# select input file
indir0 = Ldir['LOo'] + 'tef/'
# choose the tef extraction to plot
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/'

ThalMean = pickle.load(open(indir + 'ThalMean.p', 'rb'))

def plotit(ax, sect_df, sect_list, lcol):
    for sn in sect_list:
        x0 = sect_df.loc[sn,'x0']
        x1 = sect_df.loc[sn,'x1']
        y0 = sect_df.loc[sn,'y0']
        y1 = sect_df.loc[sn,'y1']
        xx = (x0+x1)/2
        yy = (y0+y1)/2
        ax.plot([x0,x1], [y0,y1], '-', color=lcol, linewidth=3)

# plotting

plt.close('all')
lw=2
fs=16
distmax = 385
fig = plt.figure(figsize=(13,8))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(223)
ax3 = fig.add_subplot(122) # map

lcol_dict = {'JdF to Strait of Georgia': 'olive',
                'JdF to Hood Canal': 'blue',
                'JdF to Whidbey Basin': 'gold',
                'JdF to South Sound': 'red'}

for lstr in lcol_dict.keys():
    sect_list, qin, qout, qsin, qsout, sin, sout, dist = ThalMean[lstr]
    lcol = lcol_dict[lstr]
    ax1.plot(dist,qin,'-o', color=lcol,linewidth=lw, label=lstr)
    #ax1.plot(dist,-qout,'-', color=lcol,linewidth=lw-1, label=lstr)
    ax1.set_xlim(0,distmax)
    ax1.grid(True)
    ax1.set_ylabel('$Q_{in}\ (1000\ m^{3}s^{-1})$', fontsize=fs)
    ax1.legend()
    if lstr == 'JdF to South Sound':
        counter = 0
        for sn in sect_list:
            sn = sn.upper()
            ax1.text(dist[counter], qin[counter]+10, sn, rotation=45, fontsize=8)
            counter += 1
        
    ax2.fill_between(dist, sout, y2=sin, color=lcol, alpha=.8)
    ax2.set_xlim(0,distmax)
    ax2.grid(True)
    ax2.set_xlabel('Distance from Mouth (km)', fontsize=fs)
    ax2.set_ylabel('$S_{in},\ S_{out}\ (psu)$', fontsize=fs)
    #ax2.legend()
    
    if lstr == 'JdF to South Sound':
        aa = [-125.5, -122, 46.7, 50.3]
        pfun.add_coast(ax3)
        pfun.dar(ax3)
        ax3.axis(aa)
        ax3.grid(True)
        
    plotit(ax3, sect_df, sect_list, lcol)

# add rivers
#%% get river info
ri_fn = Gr['ri_dir'] + 'river_info.csv'
df = pd.read_csv(ri_fn, index_col='rname')
for rn in df.index:
    try:
        fn_tr = Gr['ri_dir'] + 'tracks/' + rn + '.csv'
        df_tr = pd.read_csv(fn_tr, index_col='ind')
        x = df_tr['lon'].values
        y = df_tr['lat'].values
        ax3.plot(x, y, '-',color='purple', linewidth=2, alpha=.4)
        ax3.plot(x[-1], y[-1], '*r', alpha=.4)
        ax3.text(x[-1]+.01, y[-1]+.01, rn, alpha=.4)
    except FileNotFoundError:
        pass
        


plt.show()
    