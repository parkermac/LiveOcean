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
if False:
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
    #Litem = 'cas4_v1_lo6biom_2017.09.01_2017.09.03'
    Litem = 'cascadia1_base_lobio5_2013.01.01_2013.12.31'
print('\nProcessing ' + Litem + '\n')
Indir = indir + Litem + '/'

# USER set which section to look at
sect_name = 'shelf_45'

# set plotting parameters
if False:
    save_fig = True
    out_dir0 = Ldir['LOo'] + 'tef_plots3/'
    Lfun.make_dir(out_dir0)
    out_dir = out_dir0 + Litem + '/'
    Lfun.make_dir(out_dir, clean=True)
else:
    save_fig = False

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
    
fn = Indir + sect_name + '.nc'
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
# qq = q.mean(axis=0)
# ss = salt.mean(axis=0)
# ox = oxygen.mean(axis=0)

nday = 5
nfilt = nday*24
qq = zfun.filt_hanning_mat(q, n=nfilt)
ss = zfun.filt_hanning_mat(salt, n=nfilt)
ox = zfun.filt_hanning_mat(oxygen, n=nfilt)
# subsample
qq = qq[::nfilt, :, :]
ss = ss[::nfilt, :, :]
ox = ox[::nfilt, :, :]

ax = fig.add_subplot(221)
cs = ax.pcolormesh(xsect, z0, qq[5,:,:]/da0, vmin=-.1, vmax=.1, cmap='bwr')
fig.colorbar(cs)
ax.text(0.05, 0.1, 'Positive is ' + dir_str, transform=ax.transAxes)
ax.set_title('Mean Velocity (m/s)')

ax = fig.add_subplot(222)
cs = ax.pcolormesh(xsect, z0, ss[5,:,:], cmap='jet')
fig.colorbar(cs)
ax.set_title('Mean Salinity')

ax = fig.add_subplot(223)
cs = ax.pcolormesh(xsect, z0, ox[5,:,:], cmap='rainbow')
fig.colorbar(cs)
ax.set_title('Mean oxygen')

# add section location map
ax = fig.add_subplot(224)
plotit(ax, aa, sect_df, sect_name)


if save_fig:
    plt.savefig(out_dir + sect_name + '.png')
    plt.close()
else:
    plt.show()
