"""
Plot time mean of variables on a physical (x-or-y vs z) plot.

"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import pickle
import pandas as pd

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun
Ldir = Lfun.Lstart()

import tef_fun
from importlib import reload
reload(tef_fun)

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

indir0 = Ldir['LOo'] + 'tef/'
# choose the tef extraction to plot
run_name = Lfun.choose_item(indir0)
indir = indir0 + run_name + '/extractions/'
sect_name_nc = Lfun.choose_item(indir, tag='.nc')
fn = indir + sect_name_nc

sect_name = sect_name_nc.replace('.nc','')

# set plotting parameters
if True:
    save_fig = True
    out_dir = indir0 + run_name + '/section_plots_' + sect_name + '/'
    Lfun.make_dir(out_dir, clean=True)
else:
    save_fig = False

plt.close('all')

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

# get time axis for indexing
ot = ds['ocean_time'][:]
dt = []
for tt in ot:
    dt.append(Lfun.modtime_to_datetime(tt))
tind = np.arange(len(dt))
dt_ser = pd.Series(index=dt, data=tind)

# time variable fields
q = ds['q'][:]
salt = ds['salt'][:]

for mm in range(1,13):
    try:
        dt_mo = dt_ser[dt_ser.index.month == mm]
        it0 = dt_mo[0]
        it1 = dt_mo[-1]
    
        # form time means
        qq = q[it0:it1,:].mean(axis=0)
        ss = salt[it0:it1,:].mean(axis=0)

        fig = plt.figure(figsize=(13,8))
    
        ax = fig.add_subplot(221)
        cs = ax.pcolormesh(xsect, z0, qq/da0, vmin=-.1, vmax=.1, cmap='bwr')
        fig.colorbar(cs)
        ax.text(0.05, 0.1, 'Positive is ' + dir_str, transform=ax.transAxes)
        ax.set_title('Mean Velocity (m/s) Month = ' + str(mm))

        ax = fig.add_subplot(223)
        if 'shelf' in sect_name:
            cs = ax.pcolormesh(xsect, z0, ss, vmin = 33, vmax=34.5, cmap='rainbow')
        else:
            cs = ax.pcolormesh(xsect, z0, ss, cmap='rainbow')
        fig.colorbar(cs)
        ax.set_title('Mean Salinity')

        # add section location map
        ax = fig.add_subplot(122)
        plotit(ax, aa, sect_df, sect_name)

        if save_fig:
            nnnn = ('0000' + str(mm))[-4:]
            plt.savefig(out_dir + 'plot_' + nnnn + '.png')
            plt.close()
        else:
            plt.show()
    except IndexError:
        pass
    

