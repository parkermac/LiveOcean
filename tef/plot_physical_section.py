"""
Plot time mean of variables on a physical (x-or-y vs z) plot.

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
if True:
    print('\nSelect an Extraction to plot:\n')
    List = os.listdir(indir)
    List.sort()
    NL = len(List)
    Ldict = dict(zip(range(NL), List))
    for ii in range(NL):
        print(str(ii) + ': ' + List[ii])
    my_ii = int(input('-- Input number: '))
    Litem = Ldict[my_ii]
else:
    #Litem = 'cas4_v1_lo6biom_2017.01.01_2017.12.31'
    Litem = 'cas4_v1_lo6biom_2017.09.01_2017.09.03'
print('\nProcessing ' + Litem + '\n')
Indir = indir + Litem + '/'

LList_raw = os.listdir(indir + Litem)
LList_raw.sort()
LList = [item for item in LList_raw if '.nc' in item]

if False: # plot all .p files
    save_fig = True
    out_dir0 = Ldir['LOo'] + 'tef_plots/'
    Lfun.make_dir(out_dir0)
    out_dir = out_dir0 + Litem + '/'
    Lfun.make_dir(out_dir, clean=True)
else: # override
    save_fig = False
    LList = [item for item in LList if 'shelf_45' in item]

plt.close('all')
fig = plt.figure(figsize=(13,8))

# function for the section map
aa = [-127.4, -122, 42, 50.3]
def plotit(ax, aa, sect_df, sn):
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.grid(True)

    x0 = sect_df.loc[sn,'x0']
    x1 = sect_df.loc[sn,'x1']
    y0 = sect_df.loc[sn,'y0']
    y1 = sect_df.loc[sn,'y1']
    xx = (x0+x1)/2
    yy = (y0+y1)/2
    if (xx > aa[0]) & (xx < aa[1]) & (yy > aa[2]) & (yy < aa[3]):
        ax.plot([x0,x1], [y0,y1], '-b', linewidth=3, alpha=.5)
        ax.text(xx, yy, sn, fontsize=12, color='r', fontweight='bold')

for LL in LList:
    
    sect_name = LL.replace('.nc','')
    fn = Indir + LL
    ds = nc.Dataset(fn)
    z0 = ds['z0'][:]
    da0 = ds['DA0'][:]
    lon = ds['lon'][:]
    lat = ds['lat'][:]
    
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
        xsect = lat
    elif (x0!=x1) and (y0==y1):
        sdir = 'EW'
        if landward == 1:
            dir_str = 'Northward'
        elif landward == -1:
            dir_str = 'Southward'
        a = [x0, x1]; a.sort()
        x0 = a[0]; x1 = a[1]
        xsect = lon

    
    
    
    # time variable fields
    q = ds['q'][:]
    salt = ds['salt'][:]
    oxygen = ds['oxygen'][:]

    # form time means
    qq = q.mean(axis=0)
    ss = salt.mean(axis=0)
    ox = oxygen.mean(axis=0)
    
    ax = fig.add_subplot(221)
    cs = ax.pcolormesh(xsect, z0, qq/da0, vmin=-.2, vmax=.2, cmap='bwr')
    fig.colorbar(cs)
    ax.text(0.05, 0.1, 'Positive is ' + dir_str, transform=ax.transAxes)
    
    ax = fig.add_subplot(222)
    cs = ax.pcolormesh(xsect, z0, ss, cmap='jet')
    fig.colorbar(cs)
    
    ax = fig.add_subplot(223)
    cs = ax.pcolormesh(xsect, z0, ox, cmap='rainbow')
    fig.colorbar(cs)

    
    #ax.pcolormesh
    
    
    # add section location map
    ax = fig.add_subplot(224)
    plotit(ax, aa, sect_df, sect_name)


    if save_fig:
        plt.savefig(out_dir + sect_name + '.png')
        plt.close()
    else:
        plt.show()
