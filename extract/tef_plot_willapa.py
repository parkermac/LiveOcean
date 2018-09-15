"""
Process a TEF extraction.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import pickle

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

Indir = (Ldir['LOo'] + 'tef/' +
    'cas4_v1_lo6biom_2017.01.01_2017.12.31/')
plt.close('all')

for LL in ['willapa_mouth.p']:
    sect_name = LL.replace('.p','')
    print('\n** ' + sect_name + ' **')
    fn = Indir + LL
    Q, S, QS, qnet_lp, fnet_lp, td = tef_fun.tef_integrals(fn)

    # plotting
    lw = 2.5
    fig = plt.figure(figsize=(20,8))

    ax = plt.subplot2grid((2,3), (0,0), colspan=2)
    nt, ns = S.shape
    Td = np.tile(td.reshape(nt,1),(1, ns))
    ax.plot(td, Q[:,3]/1e3, '-r', linewidth=lw, label='Qin shallow', alpha=.5)
    ax.plot(td, Q[:,2]/1e3, '-b', linewidth=lw, label='Qout') # out
    ax.plot(td, Q[:,1]/1e3, '-r', linewidth=lw, label='Qin') # in
    ax.plot(td, Q[:,0]/1e3, '-b', linewidth=lw, label='Qout deep', alpha=.5) # out
    ax.legend(ncol=2, loc='upper left')
    ax.set_xlim(0,365)
    ax.set_ylim(-5, 5)
    ax.set_ylabel('Q (1000 m3/s)')
    ax.grid(True)
    
    ax = plt.subplot2grid((2,3), (1,0), colspan=2)
    nt, ns = S.shape
    Td = np.tile(td.reshape(nt,1),(1, ns))
    ax.plot(td, S[:,3], '-r', linewidth=lw, label='Sin shallow', alpha=.5)
    ax.plot(td, S[:,2], '-b', linewidth=lw, label='Sout') # out
    ax.plot(td, S[:,1], '-r', linewidth=lw, label='Sin') # in
    ax.plot(td, S[:,0], '-b', linewidth=lw, label='Sout deep', alpha=.5) # out
    # mark reversals
    rev = (S[:,2] - S[:,1]) > 0
    ax.plot(td[rev], 34.8*np.ones(nt)[rev], '*k', label='Reversals')
    ax.legend(loc='upper right', ncol=5)
    ax.set_xlim(0,365)
    ax.set_ylim(36, 10)
    ax.set_xlabel('Days')
    ax.set_ylabel('Salinity')
    ax.grid(True)

    # add section location map
    ax = fig.add_subplot(1,3,3)
    x0, x1, y0, y1, landward = sect_df.loc[sect_name,:]
    ax.plot([x0, x1], [y0, y1], '-m', linewidth=3)
    ax.set_title(sect_name)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-124.5, -123.6, 46, 47.4])
    
    fig.tight_layout()

    plt.show()
