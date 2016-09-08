# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 10:54:38 2016

@author: PM5

Code to check the placement of river forcing.

"""

#%% setup

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart(gridname='cascadia2', tag='frc2')

from importlib import reload
import zfun; reload(zfun)

import pfun

import netCDF4 as nc
import matplotlib.pyplot as plt

#%% where to look

f_string = 'f2013.01.01'

#%% load the river and grid

fn_riv = Ldir['LOo'] + Ldir['gtag'] + '/' + f_string + '/riv1/rivers.nc'

fn_grid = Ldir['grid'] + 'grid.nc'

dsr = nc.Dataset(fn_riv)

dsg = nc.Dataset(fn_grid)

lon_dict = dict()
lat_dict = dict()
mask_dict = dict()
tag_list = ['rho', 'u', 'v', 'psi']
for tag in tag_list:
    lon_dict[tag] = dsg.variables['lon_'+tag][:]
    lat_dict[tag] = dsg.variables['lat_'+tag][:]
    mask_dict[tag] = dsg.variables['mask_'+tag][:]

#%% plotting

plt.close()

marker_dict = {'rho': 'ok',
             'u': '>r',
             'v': '^b',
             'psi': 'xg'}

fig = plt.figure(figsize=(15,15))

ax = fig.add_subplot(111)

r_dir = dsr['river_direction'][:]
r_ix = dsr['river_Xposition'][:]
r_iy = dsr['river_Eposition'][:]
river_names = dsr['river_name'][:]
q = dsr['river_transport'][:]

for ii in range(len(r_dir)):
    rd = r_dir[ii]
    rx = r_ix[ii].astype(int)
    ry = r_iy[ii].astype(int)
    rn_arr = river_names[:,ii]
    Q = q[0, ii]

    RN = ''
    for ch in rn_arr:
        # we use .decode() to get rid of the b prefix on A
        RN = RN + ch.decode()

    mks = 25
    if rd == 0:
        if Q >= 0:
            ax.plot(lon_dict['u'][ry, rx-1], lat_dict['u'][ry, rx-1],
                    '>m', markersize=mks)
        elif Q < 0:
            ax.plot(lon_dict['u'][ry, rx-1], lat_dict['u'][ry, rx-1],
                    '<m', markersize=mks)
        ax.text(lon_dict['u'][ry, rx-1], lat_dict['u'][ry, rx-1], RN)
    elif rd == 1:
        if Q >= 0:
            ax.plot(lon_dict['v'][ry-1, rx], lat_dict['v'][ry-1, rx],
                    '^c', markersize=mks)
        elif Q < 0:
            ax.plot(lon_dict['v'][ry-1, rx], lat_dict['v'][ry-1, rx],
                    'vc', markersize=mks)
        ax.text(lon_dict['v'][ry-1, rx], lat_dict['v'][ry-1, rx], RN)

for tag in tag_list:
    ax.plot(lon_dict[tag][mask_dict[tag]==1], lat_dict[tag][mask_dict[tag]==1],
            marker_dict[tag])

pfun.add_coast(ax)
pfun.dar(ax)
ax.axis(pfun.get_aa(dsg))

plt.show()


