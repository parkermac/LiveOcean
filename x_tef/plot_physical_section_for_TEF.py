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

import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
Ldir = Lfun.Lstart()

import tef_fun

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

indir0 = Ldir['LOo'] + 'tef/'
indir = indir0 + 'cas6_v3_lo8b_2017.01.01_2017.12.31/extractions/'
sect_name = 'jdf3'
fn = indir + sect_name + '.nc'

save_fig = True

# set plotting parameters
if save_fig:
    out_dir = indir0 + 'misc_figs_cas6/'
    Lfun.make_dir(out_dir)
    out_fn = out_dir + 'TEF_physical_section_' + sect_name + '.png'
    print('Saving figure to:\n' + out_fn)

plt.close('all')
    
ds = nc.Dataset(fn)
z0 = ds['z0'][:]
da0 = ds['DA0'][:]
lon = ds['lon'][:]
lat = ds['lat'][:]

# some information about direction
x0, x1, y0, y1, landward = sect_df.loc[sect_name,:]    
if (x0==x1) and (y0!=y1):
    sdir = 'NS'
    xlab = 'Latitude'
    if landward == 1:
        dir_str = 'Eastward'
    elif landward == -1:
        dir_str = 'Westward'
    a = [y0, y1]; a.sort()
    y0 = a[0]; y1 = a[1]
    xsect = lat
elif (x0!=x1) and (y0==y1):
    sdir = 'EW'
    xlab = 'Longitude'
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

# PLOTTING
lw=3
fs=16
plt.rc('font', size=fs)

plt.close('all')

h_offset = 3

fig = plt.figure(figsize=(18,6))

for ii in [1,2,3]:

    if ii == 1:
        mm = 7
        dt_mo = dt_ser[dt_ser.index.month == mm]
        it0 = dt_mo[0]
        it1 = dt_mo[-1]
        # form time means
        qq = q[it0:it1,:,:].mean(axis=0)
        ss = salt[it0:it1,:,:].mean(axis=0)
        U = 35
    elif ii == 2:
        iit = it0+h_offset
        qq = q[iit,:,:]
        ss = salt[iit,:,:]
        U = 100
        this_dt = dt[iit]
    elif ii == 3:
        iit = it0+h_offset+6
        qq = q[iit,:,:]
        ss = salt[iit,:,:]
        U = 100
        this_dt = dt[iit]

    ax = fig.add_subplot(1,3,ii)
    cs = ax.pcolormesh(xsect, z0, 100*qq/da0, vmin=-U, vmax=U, cmap='bwr')
    ax.contour(xsect*np.ones((z0.shape[0],1)), z0, ss,
        np.arange(20,40,.2), colors='k', linewidths=.5)
    ccs = ax.contour(xsect*np.ones((z0.shape[0],1)), z0, ss,
        np.arange(20,40,1), colors='k', linewidths=2)
    ax.clabel(ccs, fmt='%d')
    
    if ii == 1:
        ax.text(0.5, 0.03, 'Positive is ' + dir_str, transform=ax.transAxes,
        ha='center', va='center', style='italic')
        ax.text(0.5, 0.08, 'Section = ' + sect_name, transform=ax.transAxes,
        ha='center', va='center', style='italic')
        ax.set_title('(a) Monthly Mean Velocity $[cm\ s^{-1}]$')
        ax.set_ylabel('Z [m]')
    elif ii in [2,3]:
        if ii == 2:
            ax.set_title('(b) Instantaneous Velocity: Flood')
        elif ii == 3:
            ax.set_title('(c) Instantaneous Velocity: Ebb')
        ax.text(.97, .03, this_dt.strftime('%Y.%m.%d %H:%M:%S'), transform=ax.transAxes,
            ha='right', va='center', style='italic')
        ax.set_yticks([])
    if sect_name == 'jdf3':
        ax.set_xticks([48.2, 48.3])
        ax.set_ylim(-200,0)
    ax.set_xlabel(xlab)
    
    
    # Inset colorbar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower left', borderpad=1) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    cb.ax.tick_params(labelsize=.8*fs)
        

fig.tight_layout()
plt.show()
plt.savefig(out_fn)
plt.rcdefaults()


