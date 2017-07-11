#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 11:06:02 2017

@author: PM5

Code to test, graphically, the results of the new make_forcing_main.py
which writes only selected surface fields to NetCDF and Azure.

"""

# setup
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
from datetime import datetime, timedelta
import netCDF4 as nc

plp = os.path.abspath('../../plotting')
if plp not in sys.path:
    sys.path.append(plp)
import pfun

import matplotlib.pyplot as plt

Ldir = Lfun.Lstart('cascadia1', 'base')
Ldir['gtagex'] = Ldir['gtag'] + '_lobio1'
f_string = 'f2017.05.18'
in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
out_name = 'ocean_surface.nc'
out_fn = in_dir + out_name

ds = nc.Dataset(out_fn)

ot = ds['ocean_time'][:]
otu = ds['ocean_time'].units

vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
vn_list2t = ['Uwind', 'Vwind', 'zeta']
vn_list3t = ['salt', 'temp', 'NO3', 'phytoplankton',
           'zooplankton', 'oxygen', 'TIC', 'alkalinity', 'PH', 'ARAG']

G = dict()
for vn in vn_list2:
    G[vn] = ds[vn][:]
    
# plotting

plt.close('all')
    
for vn in vn_list2t + vn_list3t:
    fig = plt.figure(figsize=(12,8))
    nplot = 1
    for tlev in [0, -1]:
        ax = fig.add_subplot(1,2,nplot)
        cs = ax.pcolormesh(G['lon_psi'], G['lat_psi'], ds[vn][tlev, 1:-1, 1:-1],
                           cmap='rainbow')
        try:
            tun = ds[vn].units
        except AttributeError:
            tun = ''
        ax.set_title(vn + ' (' + tun + ')')
        fig.colorbar(cs)
        pfun.dar(ax)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        t = ot[tlev]
        fs = 12
        dt = datetime(1970,1,1,0,0,0) + timedelta(days=t/86400)
        ax.text(.95, .075, dt.strftime('%Y-%m-%d'),
            horizontalalignment='right' , verticalalignment='bottom',
            transform=ax.transAxes, fontsize=fs)
        ax.text(.95, .065, dt.strftime('%H:%M') + ' UTC',
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=fs)
        nplot += 1
plt.show()
