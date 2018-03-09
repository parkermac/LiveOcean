#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This is the main program for extracting a subset of the daily output
for people in the USRS project.

Example call from the ipython command line:
    
run make_forcing_main.py -g cas3 -t v0 -x lo6m -d 2017.03.30

"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun

Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************

testing = False

from datetime import datetime
start_time = datetime.now()
import netCDF4 as nc
import zfun
import zrfun
import numpy as np

print('- Creating USRS files for ' + Ldir['date_string'])
f_string = 'f' + Ldir['date_string']

# input files
in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
fn_list_raw = os.listdir(in_dir)
fn_list = []
for item in fn_list_raw:
    if 'ocean_his' in item and '.nc' in item:
        fn_list.append(in_dir + item)
fn_list.sort()

# output location (before transfer to Azure)
out_dir00 = Ldir['LOo'] + 'usrs/'
Lfun.make_dir(out_dir00)
out_dir0 = out_dir00 + Ldir['gtagex'] + '/'
Lfun.make_dir(out_dir0)
out_dir = out_dir0 + f_string + '/'
Lfun.make_dir(out_dir, clean=True)

# get grid info
fn = fn_list[0]
ds = nc.Dataset(fn)
G, S, T = zrfun.get_basic_info(fn)
ds.close()
# get grid and indices (rho-grid) for extraction domain
x0 = -122.85
x1 = -122.6
y0 = 47.6
y1 = 47.9
Lon = G['lon_rho']
Lat = G['lat_rho']
Lon_vec = Lon[0,:]
Lat_vec = Lat[:,0]
i0, junk, junk = zfun.get_interpolant(np.array(x0), Lon_vec, extrap_nan=True)
junk, i1, junk = zfun.get_interpolant(np.array(x1), Lon_vec, extrap_nan=True)
j0, junk, junk = zfun.get_interpolant(np.array(y0), Lat_vec, extrap_nan=True)
junk, j1, junk = zfun.get_interpolant(np.array(y1), Lat_vec, extrap_nan=True)
i0 = i0[0]
i1 = i1[0] + 1
j0 = j0[0]
j1 = j1[0] + 1

lon_rho = Lon[j0:j1, i0:i1]
lat_rho = Lat[j0:j1, i0:i1]
h = G['h'][j0:j1, i0:i1]

NR, NC = lon_rho.shape
N = S['N']

fnh_list = []

if testing:
    fn_list = fn_list[:2]
    
for fn in fn_list:
    # get derived fields
    ds1 = nc.Dataset(fn)
    zeta = ds1['zeta'][0, j0:j1, i0:i1].squeeze()
    z = zrfun.get_z(h, zeta, S, only_rho=True)    
    # interpolate velocities to the rho grid
    u = (ds1['u'][0, :, j0:j1, i0-1:i1-1].squeeze()
        + ds1['u'][0, :, j0:j1, i0:i1].squeeze())/2
    v = (ds1['v'][0, :, j0-1:j1-1, i0:i1].squeeze()
        + ds1['v'][0, :, j0:j1, i0:i1].squeeze())/2   
    
    # pack into a NetCDF file:
    
    ## output files
    fnh = fn.split('/')[-1].replace('his', 'ext')
    out_fn = out_dir + fnh
    fnh_list.append(fnh)
    print(' - creating ' + fnh)
    # get rid of the old version, if it exists
    try:
        os.remove(out_fn)
    except OSError:
        pass # assume error was because the file did not exist
    ds2 = nc.Dataset(out_fn, 'w')

    vn_list2 = [ 'lon_rho', 'lat_rho', 'mask_rho', 'h']
    vn_list2p = [ 'lon_psi', 'lat_psi']
   
    vn_list2t = ['Uwind', 'Vwind', 'zeta', 'Pair',
                 'shflux', 'latent', 'sensible', 'lwrad', 'swrad']   
    vn_list3t = ['salt', 'temp', 'rho', 'w', 'AKs', 'AKv']
    
    ## Create dimensions
    ds2.createDimension('xi_rho', NC)
    ds2.createDimension('eta_rho', NR)
    ds2.createDimension('xi_psi', NC-1)
    ds2.createDimension('eta_psi', NR-1)
    ds2.createDimension('s_rho', N)
    ds2.createDimension('s_w', N+1)
    ds2.createDimension('ocean_time', None)

    # Copy variables
    vn = 'ocean_time'
    varin = ds1[vn]
    vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
    vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    vv[:] = ds1[vn][:]
        
    for vn in vn_list2:
        varin = ds1[vn]
        vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
        vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        vv[:] = ds1[vn][j0:j1, i0:i1]
    #
    for vn in vn_list2p:
        varin = ds1[vn]
        vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
        vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        vv[:] = ds1[vn][j0:j1-1, i0:i1-1]
    #
    for vn in vn_list2t:
        varin = ds1[vn]
        vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
        vv.long_name = varin.long_name
        vv.units = varin.units
        vv.time = varin.time
        vv[0,:,:] = ds1[vn][0, j0:j1, i0:i1]

    for vn in vn_list3t:
        varin = ds1[vn]
        vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
        vv.long_name = varin.long_name
        try:
            vv.units = varin.units
        except AttributeError: # salt has no units
            pass
        vv.time = varin.time
        vv[0,:,:,:] = ds1[vn][0, :, j0:j1, i0:i1]
           
    # Add derived variables
    vv = ds2.createVariable('z', float, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
    vv.long_name = 'z position'
    vv.units = 'meter'
    vv[0,:,:,:] = z
    #
    vv = ds2.createVariable('u', float, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
    vv.long_name = 'eastward velocity'
    vv.units = 'meter second-1'
    vv.time = 'ocean_time'
    vv[0,:,:,:] = u
    #
    vv = ds2.createVariable('v', float, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
    vv.long_name = 'northward velocity'
    vv.units = 'meter second-1'
    vv.time = 'ocean_time'
    vv[0,:,:,:] = v
    
    ds1.close()
    ds2.close()
    
if testing:
    # try plotting some results
    import matplotlib.pyplot as plt
    plt.close('all')
    ds = nc.Dataset(out_fn)
    plon = ds['lon_psi'][:]
    plat = ds['lat_psi'][:]
    mask = ds['mask_rho'][1:-1, 1:-1]
    pth = os.path.abspath('../../plotting/')
    if pth not in sys.path:
        sys.path.append(pth)
    import pfun
    aa = pfun.get_aa(ds)
    
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12,7))
    
    vn_list = ['salt', 'temp', 'AKs', 'u', 'v', 'w']
    
    rr = [0, 0, 0, 1, 1, 1]
    cc = [0, 1, 2, 0, 1, 2]
    
    ii = 0
    for vn in vn_list:
        ax = axes[rr[ii], cc[ii]]
        vv = ds[vn][0, -1,1:-1,1:-1].squeeze()
        vv[mask==0] = np.nan
        cs = ax.pcolormesh(plon, plat, vv, cmap='rainbow')
        fig.colorbar(cs, ax=ax)
        pfun.dar(ax)
        pfun.add_coast(ax)
        ax.set_xlim(aa[0], aa[1])
        ax.set_ylim(aa[2], aa[3])
        ax.text(.05, .95, vn, fontweight='bold', transform=ax.transAxes)
        ii += 1 
    
    ds.close()
    
    plt.show()
  
    
#%% push results to Azure
    
if not testing:

    print(' - Pushing selected files to Azure for ' + Ldir['date_string'])
    f_string = 'f' + Ldir['date_string']
    
    # Azure commands
    from azure.storage.blob import BlockBlobService
    from azure.storage.blob import PublicAccess
    ff_string = f_string.replace('.','') # azure does not like dots in container names
    # account name and key
    azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
    account = azu_dict['account']
    key = azu_dict['key']
    containername = 'usrs-' + Ldir['gtagex'].replace('_','') + '-' + ff_string
    # get a handle to the account
    blob_service = BlockBlobService(account_name=account, account_key=key)
    blob_service.create_container(containername)
    blob_service.set_container_acl(containername, public_access=PublicAccess.Container)
    
    def write_to_azure(out_fn, blob_service, containername, outname):
        # write it to Azure
        try:
            bname = open(out_fn, 'rb')
            blob_service.create_blob_from_stream(containername, out_name, bname)
            print('done putting ' + out_name)
            bname.close()
            result = 'success'
        except:
            # could be FileNotFoundError from open, or an Azure error
            print(' - Unable to write ' + out_name + ' to Azure')
            result = 'fail'
        return result
    
    result_list = []
    for out_name in fnh_list:
        out_fn = out_dir + out_name
        result_list.append(write_to_azure(out_fn, blob_service, containername, out_name))
        
    try:
        ii = result_list.index('fail')
        result = 'fail'
        print(ii)
    except ValueError: # no fails in list
        result = 'success'
      
    #%% prepare for finale
    import collections
    result_dict = collections.OrderedDict()
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
    
    #%% ************** END CASE-SPECIFIC CODE *****************
    
    ffun.finale(result_dict, Ldir, Lfun)

