"""
Calculate efflux-reflux coefficients.
"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pickle

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

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

lstr_list = ['JdF to South Sound']

for lstr in lstr_list:
    sect_list, qin, qout, qsin, qsout, sin, sout, dist = ThalMean[lstr]
    N = len(sect_list)-1
    al_up = np.nan * np.ones(N)
    al_dn = np.nan * np.ones(N)
    for ii in range(N):
        al_up[ii] = (-qout[ii]/qin[ii]) * (sout[ii]-sout[ii+1])/(sin[ii]-sout[ii+1])
        al_dn[ii] = (-qin[ii+1]/qout[ii+1]) * (sin[ii]-sin[ii+1])/(sin[ii]-sout[ii+1])
        
def plotit(ax, sect_df, sect_list, lcol):
    for sn in sect_list:
        x0 = sect_df.loc[sn,'x0']
        x1 = sect_df.loc[sn,'x1']
        y0 = sect_df.loc[sn,'y0']
        y1 = sect_df.loc[sn,'y1']
        xx = (x0+x1)/2
        yy = (y0+y1)/2
        ax.plot([x0,x1], [y0,y1], '-', color=lcol, linewidth=3)
        
def plotal(ax, sect_df, sect_list, lcol, al_up, al_dn):
    N = len(sect_list) - 1
    for ii in range(N):
        sn = sect_list[ii]
        x0 = sect_df.loc[sn,'x0']
        x1 = sect_df.loc[sn,'x1']
        y0 = sect_df.loc[sn,'y0']
        y1 = sect_df.loc[sn,'y1']
        
        SN = sect_list[ii+1]
        X0 = sect_df.loc[SN,'x0']
        X1 = sect_df.loc[SN,'x1']
        Y0 = sect_df.loc[SN,'y0']
        Y1 = sect_df.loc[SN,'y1']
        
        xx = (x0+x1 + X0+X1)/4
        yy = (y0+y1 + Y0+Y1)/4
        ax.text(xx,yy, ('%0.2f' % (al_up[ii])), color=lcol, fontsize=8, fontweight='bold')
        
plt.close('all')
lw=2
fs=16
distmax = 385
fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111) # map

lcol_dict = {'JdF to Strait of Georgia': 'olive',
                'JdF to Hood Canal': 'blue',
                'JdF to Whidbey Basin': 'gold',
                'JdF to South Sound': 'red'}

for lstr in lstr_list:
    lcol = lcol_dict[lstr]
    sect_list, qin, qout, qsin, qsout, sin, sout, dist = ThalMean[lstr]
    
    if lstr == 'JdF to South Sound':
        aa = [-125.5, -122, 47, 50.3]
        pfun.add_coast(ax)
        pfun.dar(ax)
        ax.axis(aa)
        ax.grid(True)
        
    plotit(ax, sect_df, sect_list, lcol)
    plotal(ax, sect_df, sect_list, lcol, al_up, al_dn)
    
plt.show()
    

