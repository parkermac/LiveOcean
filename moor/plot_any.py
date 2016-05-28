"""
Plots mooring records.
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
#from importlib import reload
#reload(Lfun)
import zfun

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'moor/'

# choose the type of plot to make
print('\n%s\n' % '** Choose mooring file to plot **')
m_list_raw = os.listdir(indir)
m_list = []
for m in m_list_raw:
    if '.nc' in m:
        m_list.append(m)
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
if False:
    my_npt = int(input('-- Input number -- '))
else:
    my_npt = 0 # for testing
moor_file = m_dict[my_npt]
fn = indir + moor_file

#%% calculate pH and Aragonite saturation state
import subprocess
func = ("run_co2sys(\'" +
    indir + "\',\'" +
    moor_file + "\',\'" +
    'input3' + "\')")
cmd = Ldir['which_matlab']
run_cmd = [cmd, "-nojvm", "-nodisplay", "-r", func, "&"]
proc = subprocess.Popen(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = proc.communicate() # "out" is the screen output of the matlab code
print(out.decode())

#%% load and organize data

v0_list = ['h', 'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v']
v1_list = ['ocean_time']
v2_list = []
v3_list_rho = []
v3_list_w = []
ds = nc.Dataset(fn)
for vv in ds.variables:
    vdim = ds.variables[vv].dimensions
    if ( ('ocean_time' in vdim)
        and ('s_rho' not in vdim)
        and ('s_w' not in vdim)
        and (vv != 'ocean_time') ):
        v2_list.append(vv)
    elif ( ('ocean_time' in vdim) and ('s_rho' in vdim) ):
        v3_list_rho.append(vv)
    elif ( ('ocean_time' in vdim) and ('s_w' in vdim) ):
        v3_list_w.append(vv)

# load everything into a dict
# (is this a good idea?)
V = dict()
for vv in ds.variables:
    V[vv] = ds.variables[vv][:]

ds.close()

#%% plotting

plt.close()

NR = 5; NC = 7
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(17,9), squeeze=False)

days = (V['ocean_time'] - V['ocean_time'][0])/86400.

mdays = Lfun.modtime_to_mdate_vec(V['ocean_time'])
mdt = mdates.num2date(mdays) # list of datetimes of data

# generate ticks and labels
dt_ticks = []
dt_ticks_yr = []
dt_ticklabels = []
for yr in [2013, 2014]:
    for mo in [1, 7]:
        dt = datetime(yr, mo, 1).date()
        if dt > mdt[0].date() and dt < mdt[-1].date():
            dt_ticks.append(dt)
            if mo == 1:
                dt_ticks_yr.append(dt)
            dt_ticklabels.append(dt.strftime('%Y/%m'))

cc = 0
nmid = round(V['salt'].shape[0]/2)
for vn in v3_list_rho + v3_list_w + v2_list:
    ir = int(np.floor(cc/NC))
    ic = int(cc - NC*ir)
    ax = axes[ir, ic]
    if V[vn].ndim == 2:
        ax.plot(mdt, V[vn][-1,:], '-r')
        ax.plot(mdt, V[vn][nmid,:],'-g')
        ax.plot(mdt, V[vn][8,:], '-b')
    elif V[vn].ndim == 1:
        ax.plot(mdt, V[vn])
    ax.set_xlim(mdt[0], mdt[-1])
    ax.set_xticks(dt_ticks)
    ax.set_xticklabels([])
    ax.ticklabel_format(useOffset=False, axis='y')
    ax.text(.05, .85, vn,
            horizontalalignment='left',
            transform=ax.transAxes)
    if ir == NR-1:
        ax.set_xlabel('Date')
        ax.set_xticklabels(dt_ticklabels)
    aa = ax.get_ylim()
    for dtyr in dt_ticks_yr:
        ax.plot([dtyr, dtyr], aa, '-k')
    ax.set_ylim(aa)
    cc += 1

plt.show()
