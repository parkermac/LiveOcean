#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This creates a single NetCDF file containing fields from one or more
model layers, for some time range.

NOTE: for the "vave" layer name we speed up the code by using the average (zeta=0)
layer thickness.  It takes about 1.25 minutes for a 24 hour day for the cas4 grid.

"""

from datetime import datetime, timedelta
start_time = datetime.now()
import netCDF4 as nc
import argparse

import os
import sys
pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import numpy as np
import zrfun
import zfun

# Command line arguments

# defaults
list_type = 'hourly' # hourly, daily, low_passed
layer_name = 'zeta' # name to control choices in code

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas4')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v2')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo6biom')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2017.07.20')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2017.07.22')
parser.add_argument('-lt', '--list_type', nargs='?', type=str, default=list_type)
# layer specific arguments
parser.add_argument('-ln', '--layer_name', nargs='?', type=str, default=layer_name)
args = parser.parse_args()

# save some arguments
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name
Ldir['list_type'] = args.list_type
Ldir['date_string0'] = args.date_string0
Ldir['date_string1'] = args.date_string1
Ldir['layer_name'] = args.layer_name

# make sure the output directory exists
outdir00 = Ldir['LOo']
Lfun.make_dir(outdir00)
outdir0 = outdir00 + 'layer/'
Lfun.make_dir(outdir0)
outdir = (outdir0 + Ldir['gtagex'] + '_' + Ldir['date_string0']
        + '_' + Ldir['date_string1'] + '/')
Lfun.make_dir(outdir, clean=False)

# get list of history files to plot
fn_list = Lfun.get_fn_list(args.list_type, Ldir, args.date_string0, args.date_string1)

# name output file
out_fn = (outdir + Ldir['layer_name'] + '_' + Ldir['list_type'] + '.nc')
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass
    
# lists of variables to process
dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi', 'ocean_time']
vn_list_2d = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
if args.layer_name == 'zeta':
    nlay = -1
    vn_list_2d_t = ['zeta', 'Pair']
    vn_list_3d_t = []
    vn_list_3d_t_custom = []
    vn_list_2d_uv_t = []
    vn_list_3d_uv_t = []
    vn_list_2d_custom = []
elif args.layer_name == 'bottom':
    nlay = 0
    vn_list_2d_t = ['zeta', 'Pair']
    vn_list_3d_t = ['salt', 'temp']
    vn_list_3d_t_custom = []
    vn_list_2d_uv_t = ['bustr', 'bvstr']
    vn_list_3d_uv_t = ['u', 'v']
    vn_list_2d_custom = ['zlay']
elif args.layer_name == 'svar':
    # a custom extraction of vertically integrated "mixing" in the
    # sense of destruction of salinity variance
    nlay = -1
    vn_list_2d_t = []
    vn_list_3d_t = []
    vn_list_3d_t_custom = ['mix']
    vn_list_2d_uv_t = []
    vn_list_3d_uv_t = []
    vn_list_2d_custom = ['DA']
elif args.layer_name == 'vave':
    # a custom extraction of vertically averaged properties
    nlay = -1
    vn_list_2d_t = []
    vn_list_3d_t = []
    vn_list_3d_t_custom = ['vave_salt', 'vave_temp', 'vave_rho']
    vn_list_2d_uv_t = []
    vn_list_3d_uv_t = []
    vn_list_2d_custom = []
else:
    print('Unsupported layer name')
    sys.exit()
    
# make some things
fn = fn_list[0]
G = zrfun.get_basic_info(fn, only_G=True)
DA = G['DX'] * G['DY']
ny, nx = DA.shape
h = G['h']
S = zrfun.get_basic_info(fn, only_S=True)
zr, zw = zrfun.get_z(h, 0*h, S)
zlay = zr[nlay, :, :].squeeze()

dzr = np.diff(zw, axis=0)

ds1 = nc.Dataset(fn)
ds2 = nc.Dataset(out_fn, 'w')

# Create dimensions
for dname, the_dim in ds1.dimensions.items():
    if dname in dlist:
        ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
        
# Create variables and their attributes
# - first time
vn = 'ocean_time'
varin = ds1[vn]
vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
vv.long_name = varin.long_name
vv.units = varin.units
# - then static 2d fields
for vn in vn_list_2d:
    varin = ds1[vn]
    vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
    vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    vv[:] = ds1[vn][:]
# - then static custom fields
for vn in vn_list_2d_custom:
    if vn == 'zlay':
        vv = ds2.createVariable(vn, float, ('eta_rho', 'xi_rho'))
        vv.long_name = 'Layer Z position'
        vv.units = 'm'
        vv[:] = zlay
    elif vn == 'DA':
        vv = ds2.createVariable(vn, float, ('eta_rho', 'xi_rho'))
        vv.long_name = 'Cell horizontal area'
        vv.units = 'm2'
        vv[:] = DA
#
# for time-dependent fields we create variables here and write them later
# - then time-dependent 2d fields
for vn in vn_list_2d_t:
    varin = ds1[vn]
    vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
    vv.long_name = varin.long_name
    vv.units = varin.units
    vv.time = varin.time
# - then time-dependent 3d fields
for vn in vn_list_3d_t:
    varin = ds1[vn]
    vv = ds2.createVariable(vn, varin.dtype, ('ocean_time', 'eta_rho', 'xi_rho'))
    vv.long_name = varin.long_name
    try:
        vv.units = varin.units
    except AttributeError:
        pass # salt has no units
    vv.time = varin.time
# - then time-dependent custom 3d fields (processed into 2d)
for vn in vn_list_3d_t_custom:
    if vn == 'mix':
        vv = ds2.createVariable(vn, float, ('ocean_time', 'eta_rho', 'xi_rho'))
        vv.long_name = 'Vertically Integrated Mixing'
        vv.units = 'W m-2'
        vv.time='ocean_time'
    elif vn == 'vave_salt':
        vv = ds2.createVariable(vn, float, ('ocean_time', 'eta_rho', 'xi_rho'))
        vv.long_name = 'Vertically Averaged Salinity'
        vv.time='ocean_time'
    elif vn == 'vave_temp':
        vv = ds2.createVariable(vn, float, ('ocean_time', 'eta_rho', 'xi_rho'))
        vv.long_name = 'Vertically Averaged Potential Temperature'
        vv.units = 'Celsius'
        vv.time='ocean_time'
    elif vn == 'vave_rho':
        vv = ds2.createVariable(vn, float, ('ocean_time', 'eta_rho', 'xi_rho'))
        vv.long_name = 'Vertically Averaged Density Anomaly'
        vv.units = 'kilogram meter-3'
        vv.time='ocean_time'
# - then time-dependent 2d fields interpolated from uv grids
for vn in vn_list_2d_uv_t:
    varin = ds1[vn]
    vv = ds2.createVariable(vn, varin.dtype, ('ocean_time', 'eta_rho', 'xi_rho'))
    vv.long_name = varin.long_name
    vv.units = varin.units
    vv.time = varin.time
# - then time-dependent 3d fields interpolated from uv grids
for vn in vn_list_3d_uv_t:
    varin = ds1[vn]
    vv = ds2.createVariable(vn, varin.dtype, ('ocean_time', 'eta_rho', 'xi_rho'))
    vv.long_name = varin.long_name
    vv.units = varin.units
    vv.time = varin.time
    
# Copy time dependent data
omat = np.nan * np.ones(h.shape)
omat = np.ma.masked_where(G['mask_rho']==0, omat)

tt = 0
NF = len(fn_list)
for fn in fn_list:
    if np.mod(tt,24)==0:
        print(' working on %d of %d' % (tt, NF))
        sys.stdout.flush()
    ds = nc.Dataset(fn)
    
    ds2['ocean_time'][tt] = ds['ocean_time'][0]
    
    for vn in vn_list_2d_t:
        ds2[vn][tt,:,:] = ds[vn][0, :, :]

    for vn in vn_list_3d_t:
        ds2[vn][tt,:,:] = ds[vn][0, nlay, :, :]
        
    for vn in vn_list_3d_t_custom:
        if vn == 'mix':
            zeta = ds['zeta'][0, :, :] # .squeeze() not needed unless using [:]
            salt = ds['salt'][0, :, :, :]
            zr, zw = zrfun.get_z(h, zeta, S)
            dzr = np.diff(zw, axis=0)
            # calculate net destruction of variance by vertical mixing
            K = ds['AKs'][0, 1:-1, :, :]
            dzw = np.diff(zr, axis=0)
            dsdz = np.diff(salt, axis=0) / dzw
            mix = 2 * K * dsdz**2
            #dvw = DA.reshape((1,ny,nx))*dzw
            Mix = np.sum(mix * dzw, axis=0)
            Mix = np.ma.masked_where(G['mask_rho']==0, Mix)
            ds2[vn][tt,:,:] = Mix
        elif vn == 'vave_salt':
            var = ds['salt'][0, :, :, :]
            vave_var = np.sum(var * dzr, axis=0)/h
            vave_var = np.ma.masked_where(G['mask_rho']==False, vave_var)
            ds2[vn][tt,:,:] = vave_var
        elif vn == 'vave_temp':
            var = ds['temp'][0, :, :, :]
            vave_var = np.sum(var * dzr, axis=0)/h
            vave_var = np.ma.masked_where(G['mask_rho']==False, vave_var)
            ds2[vn][tt,:,:] = vave_var
        elif vn == 'vave_rho':
            var = ds['rho'][0, :, :, :]
            vave_var = np.sum(var * dzr, axis=0)/h
            vave_var = np.ma.masked_where(G['mask_rho']==False, vave_var)
            ds2[vn][tt,:,:] = vave_var
        
    if ('sustr' in vn_list_2d_uv_t) and ('svstr' in vn_list_2d_uv_t):
        sustr0 = ds['sustr'][0, :, :]
        svstr0 = ds['svstr'][0, :, :]
        sustr = omat.copy()
        svstr = omat.copy()
        sustr[:, 1:-1] = (sustr0[:, 1:] + sustr0[:, :-1])/2
        svstr[1:-1, :] = (svstr0[1:, :] + svstr0[:-1, :])/2
        sustr = np.ma.masked_where(np.isnan(sustr), sustr)
        svstr = np.ma.masked_where(np.isnan(svstr), svstr)
        ds2['sustr'][tt,:,:] = sustr
        ds2['svstr'][tt,:,:] = svstr
        
    if ('bustr' in vn_list_2d_uv_t) and ('bvstr' in vn_list_2d_uv_t):
        bustr0 = ds['bustr'][0, :, :]
        bvstr0 = ds['bvstr'][0, :, :]
        bustr = omat.copy()
        bvstr = omat.copy()
        bustr[:, 1:-1] = (bustr0[:, 1:] + bustr0[:, :-1])/2
        bvstr[1:-1, :] = (bvstr0[1:, :] + bvstr0[:-1, :])/2
        bustr = np.ma.masked_where(np.isnan(bustr), bustr)
        bvstr = np.ma.masked_where(np.isnan(bvstr), bvstr)
        ds2['bustr'][tt,:,:] = bustr
        ds2['bvstr'][tt,:,:] = bvstr
        
    if ('u' in vn_list_3d_uv_t) and ('v' in vn_list_3d_uv_t):
        u0 = ds['u'][0, nlay, :, :]
        v0 = ds['v'][0, nlay, :, :]
        u = omat.copy()
        v = omat.copy()
        u[:, 1:-1] = (u0[:, 1:] + u0[:, :-1])/2
        v[1:-1, :] = (v0[1:, :] + v0[:-1, :])/2
        u = np.ma.masked_where(np.isnan(u), u)
        v = np.ma.masked_where(np.isnan(v), v)
        ds2['u'][tt,:,:] = u
        ds2['v'][tt,:,:] = v
        
    tt += 1
    ds.close()

ds1.close()
ds2.close()

# finale
import collections
result_dict = collections.OrderedDict()
result_dict['out_fn'] = out_fn
result_dict['date_string0'] = args.date_string0
result_dict['date_string1'] = args.date_string1
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)
if os.path.isfile(out_fn):
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'
print('')
for k in result_dict.keys():
    print('%s: %s' % (k, result_dict[k]))
    


