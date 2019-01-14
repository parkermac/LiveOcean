"""
Makes several plot that attmept to clarify the TEF
analysis method.

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
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/'
if True:
    item = Lfun.choose_item(indir, tag='.nc')
else:
    item = 'ai1.nc'
fn = indir + item

sect_name = item.replace('.nc','')

# function for the section map
aa = [-125, -122, 47, 50.3]
def plotit(ax, aa, sect_df, sn):
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.grid(True)
    x0 = sect_df.loc[sn,'x0']; x1 = sect_df.loc[sn,'x1']
    y0 = sect_df.loc[sn,'y0']; y1 = sect_df.loc[sn,'y1']
    xx = (x0+x1)/2; yy = (y0+y1)/2
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
    x_str = 'Latitude'
elif (x0!=x1) and (y0==y1):
    sdir = 'EW'
    if landward == 1:
        dir_str = 'Northward'
    elif landward == -1:
        dir_str = 'Southward'
    a = [x0, x1]; a.sort()
    x0 = a[0]; x1 = a[1]
    xsect = lon
    x_str = 'Longitude'

# time variable fields
q = ds['q'][:]
salt = ds['salt'][:]

# set time limits for plotting
daylo = 100
dayhi = 150
# daylo = 5
# dayhi = 10
hourlo = daylo*24
hourhi = dayhi*24

qq = q[hourlo:hourhi,:,:].mean(axis=0)
ss = salt[hourlo:hourhi,:,:].mean(axis=0)

plt.close('all')
fsz=14

# PLOT: Eulerian section means and Map
fig = plt.figure(figsize=(13,8))

ax = fig.add_subplot(221)
cs = ax.pcolormesh(xsect, z0, qq/da0, vmin=-.1, vmax=.1, cmap='bwr')
fig.colorbar(cs)
ax.text(0.05, 0.1, 'Positive is ' + dir_str,
    transform=ax.transAxes, fontsize=fsz)
ax.set_title('Time-Mean Velocity $(m\ s^{-1})$')
ax.set_ylabel('Z (m)')

ax = fig.add_subplot(223)
cs = ax.pcolormesh(xsect, z0, ss, cmap='jet')
Xsect = xsect.reshape((1,len(xsect))) * np.ones((ss.shape[0],1))
sci = 0.1
ax.contour(Xsect, z0, ss, np.arange(0,36,sci), colors='k', linewidths=0.5)
fig.colorbar(cs)
ax.set_title('Time-Mean Salinity (C.I. = %0.2f psu)' % (sci))
ax.set_xlabel(x_str)
ax.set_ylabel('Z (m)')

# add section location map
ax = fig.add_subplot(122)
plotit(ax, aa, sect_df, sect_name)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitue')

# PLOT: Section plots at two phases of the tide
fig = plt.figure(figsize=(13,8))

ax = fig.add_subplot(221)
cs = ax.pcolormesh(xsect, z0, q[hourlo,:,:]/da0, vmin=-1, vmax=1, cmap='bwr')
fig.colorbar(cs)
ax.text(0.05, 0.1, 'Positive is ' + dir_str, transform=ax.transAxes)
ax.set_title('Instantaneous Velocity $(m\ s^{-1})$')
ax.set_ylabel('Z (m)')

ax = fig.add_subplot(222)
cs = ax.pcolormesh(xsect, z0, q[hourlo+6,:,:]/da0, vmin=-1, vmax=1, cmap='bwr')
fig.colorbar(cs)
ax.set_title('Instantaneous Velocity $(m\ s^{-1})$ 6 hours later')

ax = fig.add_subplot(223)
cs = ax.pcolormesh(xsect, z0, salt[hourlo,:,:], cmap='jet')
ax.contour(Xsect, z0, salt[hourlo,:,:], np.arange(0,36,sci), colors='k', linewidths=0.5)
fig.colorbar(cs)
ax.set_title('Instantaneous Salinity')
ax.set_xlabel(x_str)
ax.set_ylabel('Z (m)')

ax = fig.add_subplot(224)
cs = ax.pcolormesh(xsect, z0, salt[hourlo+6,:,:], cmap='jet')
ax.contour(Xsect, z0, salt[hourlo+6,:,:], np.arange(0,36,sci), colors='k', linewidths=0.5)
fig.colorbar(cs)
ax.set_title('Instantaneous Salinity 6 hours later')
ax.set_xlabel(x_str)

    
# get detailed TEF information
# NOTE Q is packed (t,s) and s "sbinsr" going from high to zero (1000 steps)
Qi, Si, Fi, qnet_lp, fnet_lp, td, sbinsr, Q, rq, Sdiv, tef_q = tef_fun.tef_details(indir+sect_name+'.p')
smax = 36
smin = 25

# PLOT: raw and filtered transport in salinity bins
fig = plt.figure(figsize=(13,8))

ax = fig.add_subplot(211)
tef_q = np.fliplr(tef_q)
td_raw = np.arange(tef_q.shape[0])/24 # this is hourly data
tef_qmax = 40000
cs = ax.pcolormesh(td_raw[hourlo:hourhi],sbinsr,tef_q[hourlo:hourhi,:].T,
     vmin=-tef_qmax, vmax=tef_qmax, cmap='bwr')
fig.colorbar(cs, ax=ax)
ax.text(.05, .9, '(a) Hourly Transport in Salinity Bins $(m^{3}\ s^{-1}\ psu^{-1})$',
    transform=ax.transAxes, fontweight='bold', fontsize=fsz)
ax.set_ylim(smax, smin)
ax.set_ylabel('Salinity')
ax.set_xlim(daylo+2,dayhi)

ax = fig.add_subplot(212)
qmax = 4000
cs = ax.pcolormesh(td[daylo:dayhi],sbinsr,rq[daylo:dayhi,:].T, vmin=-qmax, vmax=qmax, cmap='bwr')
fig.colorbar(cs, ax=ax)
ax.text(.05, .9, '(b) $q(s)$: Tidally Averaged Transport in Salinity Bins $(m^{3}\ s^{-1}\ psu^{-1})$',
    transform=ax.transAxes, fontweight='bold', fontsize=fsz)
ax.set_ylim(smax, smin)
ax.set_ylabel('Salinity')
ax.set_xlabel('Yearday')
ax.set_xlim(daylo+2,dayhi)

# PLOT: q(s) and Q(s) with Sdiv in salinity bins

fig = plt.figure(figsize=(13,8))

ax = fig.add_subplot(211)
Qmax = Q.max()/4
cs = ax.pcolormesh(td,sbinsr,Q.T, vmin=-Qmax, vmax=Qmax, cmap='bwr')
fig.colorbar(cs, ax=ax)
ax.plot(td,Sdiv,'-k')
ax.text(.05, .9, '(a) $Q(s)$: Net Transport for Salinity >= s $(m^{3}\ s^{-1})$',
    transform=ax.transAxes, fontweight='bold', fontsize=fsz)
ax.set_ylim(smax, smin)
ax.set_ylabel('Salinity')
ax.text(.05, .1, 'Use the maximum of $Q(s)$ to find $S_{div}$ at each time',
    transform=ax.transAxes, fontweight='bold', fontsize=fsz)
ax.set_title('Black Line is the Dividing Salinity $S_{div}$', fontweight='bold', fontsize=fsz+2)

ax = fig.add_subplot(212)
qmax = rq.max()/4
cs = ax.pcolormesh(td,sbinsr,rq.T, vmin=-qmax, vmax=qmax, cmap='bwr')
fig.colorbar(cs, ax=ax)
ax.plot(td,Sdiv,'-k')
ax.text(.05, .9, '(b) $q(s)$: Tidally Averaged Transport in Salinity Bins $(m^{3}\ s^{-1}\ psu^{-1})$',
    transform=ax.transAxes, fontweight='bold', fontsize=fsz)
ax.set_ylim(smax, smin)
ax.set_ylabel('Salinity')
ax.set_xlabel('Yearday')
ax.text(.05, .1, 'Then integrate $q(s)$ above and below $S_{div}$ to find $Q_{out}$ and $Q_{in}$',
    transform=ax.transAxes, fontweight='bold', fontsize=fsz)

plt.show()
