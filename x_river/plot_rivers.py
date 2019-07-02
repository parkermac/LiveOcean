#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot as-run river time series.

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates


import os
import sys
pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import zrfun
import zfun

Ldir = Lfun.Lstart('cas6', 'v3')
fnr = Ldir['gtag'] + '_2017.01.01_2018.12.31.p'
fn = Ldir['LOo'] + 'river/' + fnr
df = pd.read_pickle(fn)

# get climatology
dfi = df.index

#cols = ['skokomish']
cols = df.columns

df_clim = pd.DataFrame(index=dfi, columns=cols)
y0 = dfi[0].year
y1 = dfi[-1].year
for riv_name in cols:
    clm_fn = Ldir['data'] + 'rivers/Data_clim/' + riv_name + '.csv'
    dfc = pd.read_csv(clm_fn, header=None, index_col=0, names=['Qr'])
    for yr in range(y0, y1+1):
        mask = dfi.year==yr
        d0 = dfi[mask][0].dayofyear
        d1 = dfi[mask][-1].dayofyear
        df_clim.loc[dfi[mask], riv_name] = dfc.loc[d0:d1,'Qr'].values

rind = df.index
dt0 = rind[0]
dt1 = rind[-1]

# plotting
plt.close('all')

for ff in range(5):
    
    fig = plt.figure(figsize=(12,10))
    
    for ii in range(9):
        rr = ii + 9*ff
        riv_name = df.iloc[:,rr].name
        ax = fig.add_subplot(3,3,ii+1)
        df.iloc[:,rr].plot(ax=ax, color='purple')
        df_clim.iloc[:,rr].plot(ax=ax, color='orange')
        ax.set_xlim(dt0, dt1)
        ax.text(.05, .9, riv_name, transform=ax.transAxes, fontweight='bold')
        if ii+1 <=6:
            ax.set_xticklabels([])
        else:
            ax.xaxis.set_major_locator(mdates.MonthLocator())
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
            ax.xaxis.set_tick_params(labelrotation=45)
            ax.set_xlabel('Date 2017')
    fig.suptitle(fnr)
    #fig.tight_layout()
        
plt.show()
