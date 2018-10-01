"""
Plots layer records.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun


Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'extract/'

# choose the mooring extraction to plot
print('\n%s\n' % '** Choose layer file to plot **')
m_list_raw = os.listdir(indir)
m_list_raw.sort()
m_list = [m for m in m_list_raw if (('.nc' in m) and ('layer_' in m))]
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
if True:
    my_npt = int(input('-- Input number -- '))
else:
    my_npt = 0 # for testing
moor_file = m_dict[my_npt]
fn = indir + moor_file

ds = nc.Dataset(fn)

# PLOTTING
plt.close('all')

fig = plt.figure(figsize=(16,8))

xp = ds['lon_psi'][:]
yp = ds['lat_psi'][:]
ot = ds['ocean_time'][:]

v_list = ['u','v', 'bustr', 'bvstr']

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


