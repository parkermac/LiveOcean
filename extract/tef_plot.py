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
print('\nProcessing ' + Litem + '\n')
LList_raw = os.listdir(indir + Litem)
LList_raw.sort()
LList = [item for item in LList_raw if ('.p' in item)]
Indir = indir + Litem + '/'

plt.close('all')

for LL in ['mb4.p']: #LList:
    
    sect_name = LL.replace('.p','')
    print('\n** ' + sect_name + ' **')

    fn = Indir + LL
    
    Q, S, QS, qnet_lp, fnet_lp, td = tef_fun.tef_integrals(fn)

    # plotting

    lw = 2.5

    fig, axes = plt.subplots(nrows=2, ncols=3, sharex=True, figsize=(20,8))

    ax = axes[0,0]
    nt, ns = S.shape
    Td = np.tile(td.reshape(nt,1),(1, ns))
    ax.plot(td, Q[:,3]/1e3, '-r', linewidth=lw, label='Qin shallow', alpha=.5)
    ax.plot(td, Q[:,2]/1e3, '-b', linewidth=lw, label='Qout') # out
    ax.plot(td, Q[:,1]/1e3, '-r', linewidth=lw, label='Qin') # in
    ax.plot(td, Q[:,0]/1e3, '-b', linewidth=lw, label='Qout deep', alpha=.5) # out
    ax.legend(ncol=2, loc='upper left')
    ax.set_xlim(0,365)
    # ax.set_ylim(-500, 500)
    #ax.set_ylim(-50, 50)
    ax.set_xlabel('Days')
    ax.set_ylabel('Q (1e3 m3/s)')
    ax.grid(True)
    

    ax = axes[1,0]
    nt, ns = S.shape
    Td = np.tile(td.reshape(nt,1),(1, ns))
    ax.plot(td, S[:,3], '-r', linewidth=lw, label='Sin shallow', alpha=.5)
    ax.plot(td, S[:,2], '-b', linewidth=lw, label='Sout') # out
    ax.plot(td, S[:,1], '-r', linewidth=lw, label='Sin') # in
    ax.plot(td, S[:,0], '-b', linewidth=lw, label='Sout deep', alpha=.5) # out
    # mark reversals
    rev = (S[:,2] - S[:,1]) > 0
    ax.plot(td[rev], 34.8*np.ones(nt)[rev], '*k', label='Reversals')
    ax.legend(loc='upper right')
    ax.set_xlim(0,365)
    #ax.set_ylim(35, 10)
    ax.set_xlabel('Days')
    ax.set_ylabel('Salinity')
    ax.grid(True)

    ax = axes[0,1]
    ax.plot(td, qnet_lp/1e3, '-k', linewidth=lw)
    ax.set_xlabel('Days')
    # ax.set_ylim(-50, 50)
    #ax.set_ylim(-5, 5)
    ax.set_xlim(0,365)
    ax.set_ylabel('LP Volume Flux (1e3 m3/s)')
    ax.grid(True)
    ax.text(0.05, 0.1, 'Positive is Landward', transform=ax.transAxes)

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
    
    x0, x1, y0, y1, landward = sect_df.loc[sect_name,:]
    ax.plot([x0, x1], [y0, y1], '-m', linewidth=3)
    ax.set_title(sect_name)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.set_xlim(-124, -122)
    ax.set_ylim(47, 49)


    plt.show()
