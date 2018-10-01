"""
Plots mooring records.

Focused on a record at the mouth of JdF.

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

# set limits
lim_dict = {'NO3': (0, 50),
        'oxygen': (0, 10),
        'u': (-1, 1)}
        
vs_dict = {'u': 'Along-channel Velocity (positive landward) $(m\ s^{-1})$',
    'NO3': 'Nitrate $(\mu M)$',
    'oxygen': 'Dissolved Oxygen $(ml\ L^{-1})$'}

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'extract/'
moor_file = 'moor_cascadia1_base_lobio5_JdFmouth_low_passed_2013.01.02_2017.12.31.nc'
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

# create a dict into which to load everything
V = dict()
# and a dict of units
Vu = dict()

#choose what to plot
#list_to_plot = v3_list_rho + v3_list_w + v2_list
list_to_plot = v3_list_rho

# hand edit variables not to look at
for v in ['CaCO3', 'PH', 'ARAG']:
    try:
        list_to_plot.remove(v)
    except ValueError:
        pass

ltp = list_to_plot.copy()
ltp.append('ocean_time')
# gather data (including ocean_time and units
for vv in ltp:
    V[vv] = ds[vv][:]
    try:
        Vu[vv] = ds[vv].units
    except AttributeError:
        Vu[vv] = ''

ds.close()

# rotate velocities to be along JdF axis
uu = V['u'].copy()
vv = V['v'].copy()
theta = -21
thrad = np.pi * theta / 180
ur = uu * np.cos(thrad) + vv * np.sin(thrad)
vr = vv * np.cos(thrad) - uu * np.sin(thrad)
V['u'] = ur
V['v'] = vr

V['oxygen'] = V['oxygen'] * 0.032/1.42903

#%% plotting

plt.close('all')

list_to_plot = ['u', 'NO3', 'oxygen']
NP = len(list_to_plot)
#NR, NC = zfun.get_rc(NP)
NR = NP
NC = 1

fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(13,8),
                         squeeze=False, sharex=True)

days = (V['ocean_time'] - V['ocean_time'][0])/86400.

mdays = Lfun.modtime_to_mdate_vec(V['ocean_time'])
mdt = mdates.num2date(mdays) # list of datetimes of data

# generate ticks and labels
dt_ticks = []
dt_ticks_yr = []
dt_ticklabels = []
yr0 = mdt[0].year
yr1 = mdt[-1].year
nyears = yr1 - yr0 + 1
ndays = (mdt[-1] - mdt[0]).days
do_ticks = False
if nyears >= 2:
    do_ticks = True
    for yr in range(yr0, yr1+1):
        for mo in [1,7]:
            dt = datetime(yr, mo, 1).date()
            if dt > mdt[0].date() and dt < mdt[-1].date():
                dt_ticks.append(dt)
                if mo == 1:
                    dt_ticks_yr.append(dt)
                if mo == 7:
                    dt_ticklabels.append(dt.strftime('%Y'))
                else:
                    dt_ticklabels.append('')
else:
    mdt = days

cc = 0
nmid = round(V['salt'].shape[0]/2)
N = V['salt'].shape[0]
nbot = 0
nmid = nmid
ntop = N-1

fs1=14
fs2=16

ir = 0
ic = 0
for vn in list_to_plot:
    
    ax = axes[ir, ic]

    nfilt = 20
    ax.plot(mdt, zfun.filt_hanning(V[vn][ntop,:],nfilt), '-r', label='Surface',linewidth=3)
    ax.plot(mdt, zfun.filt_hanning(V[vn][nmid,:],nfilt),'-g', label='Mid-Depth',linewidth=3)
    ax.plot(mdt, zfun.filt_hanning(V[vn][nbot,:],nfilt), '-b', label='Deepest',linewidth=3)
    
    ax.set_ylim(lim_dict[vn][0], lim_dict[vn][1])

    ax.tick_params(axis = 'both', which = 'major', labelsize = fs1)
    
    ax.set_xlim(mdt[0], mdt[-1])
    if do_ticks:
        ax.set_xticks(dt_ticks)
        ax.set_xticklabels([])
        if ir == NR-1:
            ax.set_xlabel('Date',fontsize=fs1)
        ax.set_xticklabels(dt_ticklabels,fontsize=fs1)
        aa = ax.get_ylim()
        dy = aa[1] - aa[0]
        for dtyr in dt_ticks_yr:
            ax.plot([dtyr, dtyr], (aa[0], aa[1]-dy/5), '-k',linewidth=3,alpha=.3)
        ax.set_ylim(aa)
    else:
        if ir == NR-1:
            ax.set_xlabel('Days')
            
    if ir==0:
        ax.set_title('LiveOcean Time Series at Mouth of Juan de Fuca',fontsize=fs2)
        ax.legend(loc='lower center',fontsize=fs2,ncol=3)

    #ax.ticklabel_format(useOffset=False, axis='y')
    fs = 16
    ax.text(.05, .85, vs_dict[vn],
            horizontalalignment='left',
            transform=ax.transAxes,
            fontsize=fs1)

            
    if ir==0:
        ax.plot([mdt[0], mdt[-1]], [0,0], '-k')
    elif ir==2:
        ax.plot([mdt[0], mdt[-1]], [1.4,1.4], '-k')
    ir += 1


plt.show()
