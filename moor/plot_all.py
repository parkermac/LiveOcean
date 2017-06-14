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
#import zfun

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
if True:
    my_npt = int(input('-- Input number -- '))
else:
    my_npt = 0 # for testing
moor_file = m_dict[my_npt]
fn = indir + moor_file

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
Vu = dict() # units

#list_to_plot = v3_list_rho + v3_list_w + v2_list
if False:
    list_to_plot = v3_list_rho
    try:
        list_to_plot.remove('CaCO3')
    except ValueError:
        pass
else: # do it by hand
    list_to_plot = ['u',
     'v',
     'temp',
     'salt',
     'NO3',
     'phytoplankton',
     'zooplankton',
     'detritus',
     'Ldetritus',
     'oxygen',
     'TIC',
     'alkalinity',
     'z_rho']#,
     # 'PH',
     # 'ARAG',
     # 'z_rho']

ltp = list_to_plot.copy()

ltp.append('ocean_time')

for vv in ltp:
    V[vv] = ds[vv][:]
    Vu[vv] = ds[vv].units
        
# an experiment about how much TIC change to expect from Ldetritus decomposition
#ad_hoc = 500
#V['dTIC'] = np.cumsum(V['Ldetritus'], axis=1)*.1*6.625 * ad_hoc
#Vu['dTIC'] = Vu['TIC']
#list_to_plot.append('dTIC')

ds.close()

#%% plotting

ylim_dict = {'alkalinity':(2200, 2540)}

NP = len(list_to_plot)

NR = np.maximum(1, np.ceil(np.sqrt(NP)).astype(int))
NC = np.ceil(np.sqrt(NP)).astype(int)

fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(13,8), squeeze=False)

days = (V['ocean_time'] - V['ocean_time'][0])/86400.

mdays = Lfun.modtime_to_mdate_vec(V['ocean_time'])
mdt = mdates.num2date(mdays) # list of datetimes of data

# generate ticks and labels
dt_ticks = []
dt_ticks_yr = []
dt_ticklabels = []
for yr in [2013, 2014, 2015, 2016]:
    for mo in [1,4,7,10]:
        dt = datetime(yr, mo, 1).date()
        if dt > mdt[0].date() and dt < mdt[-1].date():
            dt_ticks.append(dt)
            if mo == 1:
                dt_ticks_yr.append(dt)
            dt_ticklabels.append(dt.strftime('%Y/%m'))

cc = 0
#nmid = round(V['salt'].shape[0]/2)
N = V['salt'].shape[0]
nbot = 0
nmid = 3
ntop = 7
for vn in list_to_plot:
    ir = int(np.floor(cc/NC))
    ic = int(cc - NC*ir)
    ax = axes[ir, ic]
    if V[vn].ndim == 2:
        ax.plot(mdt, V[vn][ntop,:], '-r')
        ax.plot(mdt, V[vn][nmid,:],'-g')
        ax.plot(mdt, V[vn][nbot,:], '-b')
    elif V[vn].ndim == 1:
        ax.plot(mdt, V[vn])

    if True:
        # general case
        ax.set_xlim(mdt[0], mdt[-1])
        ax.set_xticks(dt_ticks)
        ax.set_xticklabels([])
        if ir == NR-1:
            ax.set_xlabel('Date')
        ax.set_xticklabels(dt_ticklabels)
        aa = ax.get_ylim()
        for dtyr in dt_ticks_yr:
            ax.plot([dtyr, dtyr], aa, '-k')
        ax.set_ylim(aa)
    else:
        # specific case
        ax.set_xlim(datetime(2015,6,1), datetime(2015,8,31))
        ax.grid()
        try:
            ax.set_ylim(ylim_dict[vn])
        except:
            pass

    ax.ticklabel_format(useOffset=False, axis='y')
    ax.text(.05, .85, vn,
            horizontalalignment='left',
            transform=ax.transAxes)
    ax.text(.05, .75, Vu[vn],
            horizontalalignment='left',
            transform=ax.transAxes)
    cc += 1

#fig.tight_layout()
plt.suptitle(fn)

plt.show()
