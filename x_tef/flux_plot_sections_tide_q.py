"""
Plot the mean of tidal energy flux and volume transport at all
TEF sections.

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
gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
Ldir = Lfun.Lstart(gridname, tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + ex_name

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

import tef_fun
import flux_fun
from importlib import reload
reload(flux_fun)

# User choices
year = 2018
year_str = str(year)
testing = True

# input
run_name = Ldir['gtagex']+'_'+year_str+'.01.01_'+year_str+'.12.31'
indir00 = Ldir['LOo'] + 'tef/'
indir0 = indir00 + run_name + '/'
indir = indir0 + 'bulk/'

# colors
clist = flux_fun.clist

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()
sect_list = list(sect_df.index)

# if testing:
#     sect_list = [sect_list[0]]
    
df = pd.DataFrame(index=sect_list)
for sn in sect_list:
    bulk = pickle.load(open(indir + sn + '.p', 'rb'))
    sx0, sx1, sy0, sy1, landward = sect_df.loc[sn,:]
    sx = (sx0+sx1)/2; sy = (sy0+sy1)/2
    if (sx0==sx1) and (sy0!=sy1):
        sdir = 'NS'
    elif (sx0!=sx1) and (sy0==sy1):
        sdir = 'EW'
    df.loc[sn,'lon'] = sx
    df.loc[sn,'lat'] = sy
    df.loc[sn,'sdir'] = sdir
    df.loc[sn,'landward'] = landward
    df.loc[sn,'F'] = bulk['fnet_lp'].mean()/1e6 # MW
    df.loc[sn,'Q'] = bulk['qnet_lp'].mean()/1e3 # 1000 m3/s
    
# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(14,8))
fs = 18
abc = 'abc'

# axis limits
x0 = -125.5; x1 = -122; y0 = 47; y1 = 50.5
dx = x1 - x0; dy = y1 - y0

ax = fig.add_subplot(121)
pfun.add_coast(ax, color='gray')
ax.axis([x0, x1, y0, y1])
pfun.dar(ax)
for sn in df.index:
    sx = df.loc[sn,'lon']
    sy = df.loc[sn,'lat']
    sdir = df.loc[sn,'sdir']
    landward = df.loc[sn,'landward']
    F = df.loc[sn,'F']
    scl = 3e4
    if sdir == 'EW':
        ax.quiver(sx,sy, 0, landward*F, scale=scl, scale_units='height')
    elif sdir == 'NS':
        ax.quiver(sx,sy, landward*F, 0, scale=scl, scale_units='height')

ax = fig.add_subplot(122)
pfun.add_coast(ax, color='gray')
ax.axis([x0, x1, y0, y1])
pfun.dar(ax)
ax.set_yticklabels([])
for sn in df.index:
    sx = df.loc[sn,'lon']
    sy = df.loc[sn,'lat']
    sdir = df.loc[sn,'sdir']
    landward = df.loc[sn,'landward']
    Q = df.loc[sn,'Q']
    scl = 100
    if sdir == 'EW':
        ax.quiver(sx,sy, 0, landward*Q, scale=scl, scale_units='height')
    elif sdir == 'NS':
        ax.quiver(sx,sy, landward*Q, 0, scale=scl, scale_units='height')
        
plt.show()


    

    