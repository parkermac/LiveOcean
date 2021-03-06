# -*- coding: utf-8 -*-
"""
Extract a CTD-like cast at a single time and place

"""

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import numpy as np
from datetime import datetime
import zfun
import zrfun
import netCDF4 as nc

# set defaults
gridname = 'cas3'
tag = 'v0'
ex_name = 'lo6m'
date_string = datetime(2017,3,31).strftime(format='%Y.%m.%d')
sta_name = 'HCB010'
lon_str = '-122.8200'
lat_str = '47.6670'

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# optional input arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default=gridname)
parser.add_argument('-t', '--tag', nargs='?', type=str, default=tag)
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default=ex_name)
parser.add_argument('-s', '--sta_name', nargs='?', type=str, default=sta_name)
parser.add_argument('-d', '--date_string', nargs='?', type=str, default=date_string)
parser.add_argument('-lon', '--lon_str', nargs='?', type=str, default=lon_str)
parser.add_argument('-lat', '--lat_str', nargs='?', type=str, default=lat_str)
args = parser.parse_args()

# get the dict Ldir
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name
Ldir['date_string'] = args.date_string
Ldir['sta_name'] = args.sta_name
Ldir['lon_str'] = args.lon_str
Ldir['lat_str'] = args.lat_str

# make sure the output directory exists
outdir = Ldir['LOo'] + 'cast/'
Lfun.make_dir(outdir)

dt = Ldir['date_string']

#%% function definitions

def get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy):
    dims = ds.variables[vv].dimensions
    if 'eta_rho' in dims:
        grd = 'rho'
    elif 'eta_u' in dims:
        grd = 'u'
    elif 'eta_v' in dims:
        grd = 'v'
    else:
        print('grid error!')
    xi0 = Xi0[grd]; yi0 = Yi0[grd]
    xi1 = Xi1[grd]; yi1 = Yi1[grd]
    aix = Aix[grd]; aiy = Aiy[grd]

    xi01 = np.array([xi0, xi1]).flatten()
    yi01 = np.array([yi0, yi1]).flatten()
    return xi01, yi01, aix, aiy

#%% set up for the extraction


# target position
Lon = np.array(float(Ldir['lon_str']))
Lat = np.array(float(Ldir['lat_str']))

# get grid info
indir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + dt + '/'
fn = indir + 'ocean_his_0021.nc' # noon local standart time
G = zrfun.get_basic_info(fn, only_G=True)
S = zrfun.get_basic_info(fn, only_S=True)

# get interpolants for this point
Xi0 = dict(); Yi0 = dict()
Xi1 = dict(); Yi1 = dict()
Aix = dict(); Aiy = dict()
for grd in ['rho', 'u', 'v']:
    xx = G['lon_' + grd][1,:]
    yy = G['lat_' + grd][:,1]
    xi0, xi1, xfr = zfun.get_interpolant(Lon, xx, extrap_nan=True)
    yi0, yi1, yfr = zfun.get_interpolant(Lat, yy, extrap_nan=True)
    Xi0[grd] = xi0
    Yi0[grd] = yi0
    Xi1[grd] = xi1
    Yi1[grd] = yi1
    # create little arrays that are used in the actual interpolation
    Aix[grd] = np.array([1-xfr, xfr]).reshape((1,1,2))
    Aiy[grd] = np.array([1-yfr, yfr]).reshape((1,2))

# generating some lists
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

V = dict()
V_long_name = dict()
V_units = dict()
v_all_list = v0_list + v1_list + v2_list + v3_list_rho + v3_list_w
for vv in v_all_list:
    V[vv] = np.array([])
    try:
        V_long_name[vv] = ds.variables[vv].long_name
    except:
        V_long_name[vv] = ''
    try:
        V_units[vv] = ds.variables[vv].units
    except:
        V_units[vv] = ''

# get static variables
for vv in v0_list:
    xi01, yi01, aix, aiy = get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
    vvtemp = ds.variables[vv][yi01, xi01].squeeze()
    V[vv] =   ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)

ds.close()

#%% extract time-dependent fields


print('Working on date_list item: ' + dt)
sys.stdout.flush()

ds = nc.Dataset(fn)
for vv in v1_list:
    vtemp = ds.variables[vv][:].squeeze()
    V[vv] = np.append(V[vv], vtemp)
for vv in v2_list:
    xi01, yi01, aix, aiy = get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
    vvtemp = ds.variables[vv][:, yi01, xi01].squeeze()
    vtemp = ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
    V[vv] = np.append(V[vv], vtemp)
for vv in v3_list_rho:
    xi01, yi01, aix, aiy = get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
    vvtemp = ds.variables[vv][:, :, yi01, xi01].squeeze()
    vtemp = ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
    V[vv] = vtemp.reshape((S['N'],1))
for vv in v3_list_w:
    xi01, yi01, aix, aiy = get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
    vvtemp = ds.variables[vv][:, :, yi01, xi01].squeeze()
    vtemp = ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
    V[vv] = vtemp.reshape((S['N']+1,1))

ds.close()

# create z_rho and z_w (has to be done after we have V['zeta'])
hh = V['h'][:] * np.ones_like(V['zeta'])
z_rho, z_w = zrfun.get_z(hh, V['zeta'][:], S)
V['hh'] = hh
V_long_name['hh'] = 'bottom depth (positive down) as a vector'
V_units['hh'] = 'm'
V['z_rho'] = z_rho
V_long_name['z_rho'] = 'z on rho points (positive up)'
V_units['z_rho'] = 'm'
V['z_w'] = z_w
V_long_name['z_w'] = 'z on w points (positive up)'
V_units['z_w'] = 'm'
v2_list.append('hh')
v3_list_rho.append('z_rho')
v3_list_w.append('z_w')

#%% save the output to NetCDF
out_fn = (outdir +
    Ldir['gtagex'] + '_' +
    Ldir['sta_name'] + '_' +
    Ldir['date_string'] +
    '.nc')
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
foo = nc.Dataset(out_fn, 'w')

N = S['N']
NT = len(V['ocean_time'][:])

foo.createDimension('scalar', 1)
foo.createDimension('s_rho', N)
foo.createDimension('s_w', N+1)
foo.createDimension('ocean_time', NT)

for vv in v0_list:
    v_var = foo.createVariable(vv, float, ('scalar'))
    v_var[:] = V[vv][:]
    v_var.long_name = V_long_name[vv]
    v_var.units = V_units[vv]
for vv in v1_list:
    v_var = foo.createVariable(vv, float, ('ocean_time',))
    v_var[:] = V[vv][:]
    v_var.long_name = V_long_name[vv]
    v_var.units = V_units[vv]
for vv in v2_list:
    v_var = foo.createVariable(vv, float, ('ocean_time',))
    v_var[:] = V[vv][:]
    v_var.long_name = V_long_name[vv]
    v_var.units = V_units[vv]
for vv in v3_list_rho:
    v_var = foo.createVariable(vv, float, ('s_rho', 'ocean_time'))
    v_var[:] = V[vv][:]
    v_var.long_name = V_long_name[vv]
    v_var.units = V_units[vv]
for vv in v3_list_w:
    v_var = foo.createVariable(vv, float, ('s_w', 'ocean_time'))
    v_var[:] = V[vv][:]
    v_var.long_name = V_long_name[vv]
    v_var.units = V_units[vv]
    
# testing
    
salt = foo['salt'][:].squeeze()
temp = foo['temp'][:].squeeze()
z = foo['z_rho'][:].squeeze()

import matplotlib.pyplot as plt
fig = plt.figure(figsize=(12,7))

ax = fig.add_subplot(121)
ax.plot(salt, z)

ax = fig.add_subplot(122)
ax.plot(temp, z)

plt.show()

foo.close()

