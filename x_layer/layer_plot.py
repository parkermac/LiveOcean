"""
Plots layer records.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime

import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun


Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'layer/'

# choose the layer extraction to plot
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/'
infile = Lfun.choose_item(indir, tag='.nc')
fn = indir + infile

ds = nc.Dataset(fn)

# PLOTTING
plt.close('all')

fig = plt.figure(figsize=(16,8))

xp = ds['lon_psi'][:]
yp = ds['lat_psi'][:]
ot = ds['ocean_time'][:]

v_list = ['u','v', 'bustr', 'bvstr', 'ubar', 'vbar']
#v_list = ['zeta', 'Pair']
#v_list = ['vave_salt', 'vave_temp', 'vave_rho']
#v_list = ['salt', 'temp', 'oxygen']

F = dict() # field
S = dict() # series
for vn in v_list:
    F[vn] = ds[vn][0,:,:].squeeze()
    S[vn] = ds[vn][:,10,10].squeeze()

days = (ot - ot[0])/86400.
mdays = Lfun.modtime_to_mdate_vec(ot)
mdt = mdates.num2date(mdays) # list of datetimes of data

NC = len(v_list)
count = 1
for vn in v_list:
    ax = fig.add_subplot(2,NC,count)
    cs = ax.pcolormesh(xp, yp, F[vn][1:-1, 1:-1])
    fig.colorbar(cs, ax=ax)
    ax.set_title(vn)
    count += 1

for vn in v_list:
    ax = fig.add_subplot(2,NC,count)
    ax.plot(days, S[vn])
    ax.set_title(vn)
    count += 1


plt.suptitle(fn)

plt.show()


