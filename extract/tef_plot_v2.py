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

indir = Ldir['LOo'] + 'tef/'
if False:
    print('\nSelect an Extraction to plot:\n')
    List_raw = os.listdir(indir)
    List_raw.sort()
    List = [item for item in List_raw]
    NL = len(List)
    Ldict = dict(zip(range(NL), List))
    for ii in range(NL):
        print(str(ii) + ': ' + List[ii])
    if True:
        my_ii = int(input('-- Input number: '))
    else:
        my_ii = 0 # for testing
    Litem = Ldict[my_ii]
else:
    Litem = 'cas4_v1_lo6biom_2017.01.01_2017.12.31'
    
print('\nProcessing ' + Litem + '\n')
LList_raw = os.listdir(indir + Litem)
LList_raw.sort()
LList = [item for item in LList_raw if ('.p' in item) and ('sog' in item)]
Indir = indir + Litem + '/'

plt.close('all')

for LL in LList:
    
    sect_name = LL.replace('.p','')
    print('\n** ' + sect_name + ' **')

    fn = Indir + LL
    
    Qi, Si, Fi, qnet_lp, fnet_lp, td = tef_fun.tef_integrals_v2(fn)

    # plotting

    lw = 2.5

    fig, axes = plt.subplots(nrows=2, ncols=3, sharex=True, figsize=(20,8))

    ax = axes[0,0]
    nt, ns = Si.shape
    Td = np.tile(td.reshape(nt,1),(1, ns))
    ax.plot(td, Qi[:,1]/1e3, '-b', linewidth=lw, label='Q1') 
    ax.plot(td, Qi[:,0]/1e3, '-r', linewidth=lw, label='Q0') 
    ax.legend(ncol=2, loc='upper left')
    ax.set_xlim(0,365)
    # ax.set_ylim(-500, 500)
    #ax.set_ylim(-50, 50)
    ax.set_xlabel('Days')
    ax.set_ylabel('Q (1e3 m3/s)')
    ax.grid(True)
    

    ax = axes[1,0]
    Td = np.tile(td.reshape(nt,1),(1, ns))
    ax.plot(td, Si[:,1], '-b', linewidth=lw, label='S1')
    ax.plot(td, Si[:,0], '-r', linewidth=lw, label='S0')
    # mark reversals
    rev = (Qi[:,1] - Qi[:,0]) > 0
    ax.plot(td[rev], 34.8*np.ones(nt)[rev], '*k', label='Reversals')
    ax.legend(loc='lower right', ncol=2)
    ax.set_xlim(0,365)
    #ax.set_ylim(35, 10)
    ax.set_xlabel('Days')
    ax.set_ylabel('Salinity')
    ax.grid(True)

    ax = axes[0,1]
    ax.plot(td, qnet_lp/1e3, '-k', linewidth=lw)
    ax.plot(td, np.nansum(Qi, axis=1)/1e3, '-m', linewidth=lw-1, alpha=.5)
    ax.set_xlabel('Days')
    # ax.set_ylim(-50, 50)
    #ax.set_ylim(-5, 5)
    ax.set_xlim(0,365)
    ax.set_ylabel('LP Volume Flux (1e3 m3/s)')
    ax.grid(True)
    
    # some information about direction
    x0, x1, y0, y1, landward = sect_df.loc[sect_name,:]    
    if (x0==x1) and (y0!=y1):
        sdir = 'NS'
        if landward == 1:
            dir_str = 'Eastward'
        elif landward == -1:
            dir_str = 'Westward'
        a = [y0, y1]; a.sort()
        y0 = a[0]; y1 = a[1]
    elif (x0!=x1) and (y0==y1):
        sdir = 'EW'
        if landward == 1:
            dir_str = 'Northward'
        elif landward == -1:
            dir_str = 'Southward'
        a = [x0, x1]; a.sort()
        x0 = a[0]; x1 = a[1]
    ax.text(0.05, 0.1, 'Positive is ' + dir_str, transform=ax.transAxes)

    ax = axes[1,1]
    ax.plot(td, fnet_lp/1e9, '-k', linewidth=lw)
    ax.set_xlim(0,365)
    ax.set_xlabel('Days')
    #ax.set_ylim(-15, 15)
    #ax.set_ylim(-1, 1)
    ax.set_ylabel('Energy Flux (GW)')
    ax.grid(True)

    # remove extra axes
    axes[0,2].set_axis_off()
    axes[1,2].set_axis_off()

    # add section location map
    ax = fig.add_subplot(1,3,3)
    
    ax.plot([x0, x1], [y0, y1], '-m', linewidth=3)
    ax.set_title(sect_name)
    pfun.add_coast(ax)
    pfun.dar(ax)
    if 'willapa' in sect_name:
        ax.set_xlim(-124.5, -123.6)
        ax.set_ylim(46, 47.4)
    else:
        ax.set_xlim(-124, -122)
        ax.set_ylim(47, 49)


    plt.show()
