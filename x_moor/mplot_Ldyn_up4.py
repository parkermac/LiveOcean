"""
Plot a mooring record just focused on physical variables for the
tracker experiments in EJdF3d_3d_up4.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import pandas as pd

import os, sys
sys.path.append(os.path.abspath('../alpha'))

import Lfun
import zfun

Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'moor/'

if False:
    # choose the mooring extraction to plot
    item = Lfun.choose_item(indir0)
    indir = indir0 + item + '/'
    infile = Lfun.choose_item(indir, tag='.nc')
    fn = indir + infile
else:
    fn = indir0 + 'cas6_v3_lo8da_2018.05.15_2018.06.15/AISouth_hourly.nc'
ds = nc.Dataset(fn)

# select indices of depths to plot
N = ds['z_rho'][:].shape[1]
nlist = [N-1, 0]
# and associated colors
cdict = dict(zip(nlist, ['orange','cornflowerblue']))

# get fields
ubar = ds['ubar'][:]
vbar = ds['vbar'][:]
u = ds['u'][:]
v = ds['v'][:]
salt = ds['salt'][:]
temp = ds['temp'][:]
z_rho = ds['z_rho'][:]
ot = ds['ocean_time'][:]
t_md = Lfun.modtime_to_mdate_vec(ot)
t_dt = []
for t in ot:
    t_dt.append(Lfun.modtime_to_datetime(t))
    
def rot_vec(u,v,theta):
    ur = u*np.cos(theta) + v*np.sin(theta)
    vr = v*np.cos(theta) - u*np.sin(theta)
    return ur, vr
    
def de_mean(u, v):
    up = u - np.nanmean(u)
    vp = v - np.nanmean(v)
    return up, vp
    
# find angle of principal axes
up, vp = de_mean(ubar, vbar)
theta = 0.5 * np.arctan2(2*np.nanmean(up*vp),(np.nanvar(up)-np.nanvar(vp)))

# and rotate
ubar_r, vbar_r = rot_vec(ubar, vbar, theta)
u_r, v_r = rot_vec(u,v,theta)

# PLOTTING
plt.close('all')
fs = 16
lw = 2
plt.rc('font', size=fs)
fig = plt.figure(figsize=(16,12))

ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)
for nn in nlist:
    z = int(z_rho[:,nn].mean())
    df = pd.DataFrame(index=t_dt,columns=['u','v','ulp','vlp','salt','temp'])
    df['u'] = u_r[:,nn]
    df['v'] = v_r[:,nn]
    df['ulp'] = zfun.filt_godin(u_r[:,nn])
    df['vlp'] = zfun.filt_godin(v_r[:,nn])
    df['salt'] = salt[:,nn]
    df['temp'] = temp[:,nn]
    df['u'].plot(ax=ax1, ls='-', color=cdict[nn], label='Z = '+str(z)+' m',legend=False, grid=True)
    df['ulp'].plot(ax=ax1, ls='--', color=cdict[nn], label='Z = '+str(z)+' m',legend=False, grid=True)
    if nn == 0:
        # add boxes for time periods
        cclist = ['r', 'b', 'r', 'b']
        ii = 0
        for rds in ['2018.05.15', '2018.05.22', '2018.05.29', '2018.06.05']:
            rdt0 = datetime.strptime(rds, '%Y.%m.%d')
            rdt1 = rdt0 + timedelta(days=8)
            df['u'][rdt0:rdt1].plot(ax=ax1, ls='-', lw=5, color=cclist[ii],alpha=.5)
            ii += 1
    
    df['salt'].plot(ax=ax2, ls='-', color=cdict[nn], label='Z = '+str(z)+' m',legend=True, grid=True)
    df['temp'].plot(ax=ax3, ls='-', color=cdict[nn], label='Z = '+str(z)+' m',legend=False, grid=True)
    ax1.set_xlim(t_dt[0], t_dt[-1])
    ax2.set_xlim(t_dt[0], t_dt[-1])
    ax3.set_xlim(t_dt[0], t_dt[-1])
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax1.set_xticklabels([], minor=True)
    ax2.set_xticklabels([], minor=True)
    
    ax1.text(.05,.9,r'Along-channel Velocity $[m\ s^{-1}]$', transform=ax1.transAxes)
    ax2.text(.05,.9,r'Salinity $[g\ kg^{-1}]$', transform=ax2.transAxes)
    ax3.text(.05,.9,r'Potential Temperature $[^{\circ}C]$', transform=ax3.transAxes)
    


plt.show()
plt.rcdefaults()


