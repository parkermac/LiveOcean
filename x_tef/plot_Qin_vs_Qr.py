"""
Compare Qin/out to Qr at a chosen location.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
import pickle
import pandas as pd

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun

import tef_fun
from importlib import reload
reload(tef_fun)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

Ldir = Lfun.Lstart()

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

indir = Ldir['LOo'] + 'tef/'
Litem = 'cas4_v1_lo6biom_2017.01.01_2017.12.31'
print('\nProcessing ' + Litem + '\n')
Indir = indir + Litem + '/'

save_fig = True
out_dir = Indir + 'plots/'
Lfun.make_dir(out_dir)


sr_list = [('sog1', 'fraser'),
            ('jdf2', 'fraser'),
            ('ai3', 'skagit')]


for sr in sr_list:

    sect_name = sr[0]
    riv_name = sr[1]

    # get some river data
    riv_indir = Ldir['parent'] + 'ptools_output/river/all_2017/'
    riv_fn = riv_indir + riv_name + '.p'
    riv = pickle.load(open(riv_fn, 'rb'))

    # get TEF records
    tef_fn = Indir + sect_name + '.p'
    Qi, Si, Fi, qnet_lp, fnet_lp, td = tef_fun.tef_integrals(tef_fn)

    # convert time of the riv Series to day 
    tdr = np.array(riv.index.dayofyear) + 0.5
    Qr = riv.values

    df = pd.DataFrame(index=(np.arange(366)+.5))
    tdi = pd.Index(td)
    df.loc[tdi,'Qin'] = Qi[:,0]/1e3
    df.loc[tdi,'Qout'] = -Qi[:,1]/1e3
    df.loc[tdi,'Sin'] = Si[:,0]
    df.loc[tdi,'Sout'] = Si[:,1]
    df.loc[tdi, 'DS'] = df.Sin-df.Sout
    tdri = pd.Index(tdr)
    df.loc[tdri,'Qr'] = Qr/1e3

    # plotting
    fs = 14
    fig =  plt.figure(figsize=(16,10))

    ax = fig.add_subplot(2,3,1)
    df.plot(y=['Qin'], ax=ax, legend=False)
    ax.set_xlim(0,365)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Yearday 2017')
    ax.set_ylabel('Qin [1000 m3/s]')
    ax.grid(True)
    ax.text(.05, .9, 'Qin',
        horizontalalignment='left',
        transform=ax.transAxes, fontweight='bold', fontsize=fs)

    ax = fig.add_subplot(2,3,4)
    df.plot(y=['Sin','Sout'], ax=ax, legend=False)
    ax.set_xlim(0,365)
    ax.grid(True)
    ax.set_xlabel('Yearday 2017')
    ax.set_ylabel('Salinity (psu)')
    ax.text(.05, .9, 'Sin and Sout',
        horizontalalignment='left',
        transform=ax.transAxes, fontweight='bold', fontsize=fs)

    ax = fig.add_subplot(2,3,2)
    dfn = df.dropna()
    df.plot(x='Qr', y='Qin', ax=ax, style='*', legend=False)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Qr [1000 m3/s]')
    ax.set_ylabel('Qin [1000 m3/s]')
    ax.grid(True)
    ax.text(.05, .9, 'Qin vs. Qr',
        horizontalalignment='left',
        transform=ax.transAxes, fontweight='bold', fontsize=fs)

    ax = fig.add_subplot(2,3,5)
    dfn = df.dropna()
    df.plot(x='Qr', y='DS', ax=ax, style='*', legend=False)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Qr [1000 m3/s]')
    ax.set_ylabel('Salinity (psu)')
    ax.grid(True)
    ax.text(.05, .9, '(Sin - Sout) vs. Qr',
        horizontalalignment='left',
        transform=ax.transAxes, fontweight='bold', fontsize=fs)

    ax = fig.add_subplot(2,3,3)
    df.plot(y='Qr', ax=ax, legend=False)
    ax.set_xlim(0,365)
    ax.set_xlabel('Yearday 2017')
    ax.set_ylabel('Qr [1000 m3/s]')
    ax.grid(True)
    ax.text(.05, .9, riv_name.capitalize(),
        horizontalalignment='left',
        transform=ax.transAxes, fontweight='bold', fontsize=fs)

    # add section location map
    x0, x1, y0, y1, landward = sect_df.loc[sect_name,:]    
    ax = fig.add_subplot(2,3,6)
    ax.plot([x0, x1], [y0, y1], '-m', linewidth=3)
    pfun.add_coast(ax)
    pfun.dar(ax)
    aa = [-125.5, -122, 47, 49.5]
    ax.axis(aa)
    ax.text(.05, .9, 'Transport section: ' + sect_name,
        horizontalalignment='left',
        transform=ax.transAxes, fontweight='bold', fontsize=fs)

    if save_fig:
        out_fn = 'TEF_Qr_' + sect_name + '_' + riv_name + '.png'
        print(out_fn)
        plt.savefig(out_dir + out_fn)
        plt.close()
    else:
        plt.show()
