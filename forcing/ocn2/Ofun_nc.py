#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 13:51:06 2017

@author: PM5

Functions for writing saved fields to ROMS NetCDF.
"""

import os
import netCDF4 as nc
import numpy as np
import Lfun

def make_clm_file(Ldir, nc_dir, fh_dir, c_dict, dt_list, S, G):
    # name output file
    clm_fn = nc_dir + 'ocean_clm.nc'
    # get rid of the old version, if it exists
    try:
        os.remove(clm_fn)
    except OSError:
        pass # assume error was because the file did not exist
    foo = nc.Dataset(clm_fn, 'w', format='NETCDF3_CLASSIC')
    # create dimensions
    for vn in ['salt', 'temp', 'v3d', 'v2d', 'zeta', 'ocean']:
        foo.createDimension(vn+'_time', len(dt_list)) # could use None
    foo.createDimension('s_rho', S['N'])
    for tag in ['rho', 'u', 'v']:
        foo.createDimension('eta_'+tag, G['lat_'+tag].shape[0])
        foo.createDimension('xi_'+tag, G['lon_'+tag].shape[1])
    # add time data
    dtm_list = []
    for dt in dt_list:
        dtm_list.append(Lfun.datetime_to_modtime(dt))        
    for vn in ['salt', 'temp', 'v3d', 'v2d', 'zeta', 'ocean']:
        vv = foo.createVariable(vn+'_time', float, (vn+'_time',))
        vv.units = 'seconds since 1970.01.01 UTC'
        vv[:] = np.array(dtm_list)
    # add 2d field data
    vv = foo.createVariable('zeta', float, ('zeta_time', 'eta_rho', 'xi_rho'))
    vv.long_name = 'sea surface height climatology'
    vv.units = 'meter'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:] = c['ssh']
    vv = foo.createVariable('ubar', float, ('v2d_time', 'eta_u', 'xi_u'))
    vv.long_name = 'vertically averaged u-momentum climatology'
    vv.units = 'meter second-1'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:] = c['ubar']
    vv = foo.createVariable('vbar', float, ('v2d_time', 'eta_v', 'xi_v'))
    vv.long_name = 'vertically averaged v-momentum climatology'
    vv.units = 'meter second-1'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:] = c['vbar']
    # add 3d field data
    vv = foo.createVariable('u', float, ('v3d_time', 's_rho', 'eta_u', 'xi_u'))
    vv.long_name = 'u-momentum component climatology'
    vv.units = 'meter second-1'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:,:] = c['u3d']
    vv = foo.createVariable('v', float, ('v3d_time', 's_rho', 'eta_v', 'xi_v'))
    vv.long_name = 'v-momentum component climatology'
    vv.units = 'meter second-1'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:,:] = c['v3d']
    vv = foo.createVariable('salt', float, ('salt_time', 's_rho', 'eta_rho', 'xi_rho'))
    vv.long_name = 'salinity climatology'
    vv.units = 'PSU'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:,:] = c['s3d']
    vv = foo.createVariable('temp', float, ('temp_time', 's_rho', 'eta_rho', 'xi_rho'))
    vv.long_name = 'potential temperature climatology'
    vv.units = 'Celsius'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:,:] = c['theta']
    foo.close()
    print('-Writing ocean_clm.nc')

def add_bio(nc_dir):
    # name output file
    clm_fn = nc_dir + 'ocean_clm.nc'
    foo = nc.Dataset(clm_fn, 'a', format='NETCDF3_CLASSIC')
#    vn_dict =  {'NO3':'MicroMolar',
#                'phytoplankton':'MicroMolar N',
#                'zooplankton':'MicroMolar N',
#                'detritus':'MicroMolar N',
#                'Ldetritus':'MicroMolar N',
#                'CaCO3':'MicroMolar C',
#                'oxygen':'MicroMolar O',
#                'alkalinity':'MicroMolar',
#                'TIC':'MicroMolar C'}
    vn_dict =  {'NO3':'MicroMolar'}

    salt = foo['salt'][:]
    for vn in vn_dict.keys():
        vv = foo.createVariable(vn, float, ('salt_time', 's_rho', 'eta_rho', 'xi_rho'))
        vv.long_name = vn + ' climatology'
        vv.units = vn_dict[vn]
        vv.time = 'salt_time'
        V = create_bio_var(salt, vn)
        vv[:] = V
    foo.close()
    print('-Writing bio variables to ocean_clm.nc')

def create_bio_var(salt, vn):
    if vn == 'NO3':
        # Salinity vs. NO3 [uM], Ryan McCabe 8/2015
        # NO3 = mm*salt + bb
        mm = 0*salt
        bb = 0*salt
        ind = (salt < 31.898)
        mm[ind] = 0
        bb[ind] = 0
        ind = ((salt >= 31.898) & (salt < 33.791))
        mm[ind] = 16.3958
        bb[ind] = -522.989
        ind = ((salt >= 33.791) & (salt < 34.202))
        mm[ind] = 29.6973
        bb[ind] = -972.4545
        ind = ((salt >= 34.202) & (salt < 34.482))
        mm[ind] = 8.0773
        bb[ind] = -233.0007
        ind = ((salt >= 34.482) & (salt < 35))
        mm[ind] = -28.6251
        bb[ind] = 1032.5686
        NO3 = mm*salt + bb
        # Set maximum NO3 to 45 microMolar (found at ~800m depth), based on
        # evidence from historical NO3 data in NODC World Ocean Database.
        NO3[NO3 > 45] = 45
        # Ensure that there are no negative values.
        NO3[NO3 < 0] = 0
        return NO3
    elif vn == ['oxygen']:
        if np.nanmax(salt) > 36:
            print('Salt out of range for oxgen regression')
        # Salinity vs. oxygen [uM], Ryan McCabe 8/2015
        # oxygen = mm*salt + bb
        mm = 0*salt
        bb = 0*salt
        ind = (salt < 32.167)
        mm[ind] = 0
        bb[ind] = 300
        ind = ((salt >= 32.167) & (salt < 33.849))
        mm[ind] = -113.9481
        bb[ind] = 3965.3897
        ind = ((salt >= 33.849) & (salt < 34.131))
        mm[ind] = -278.3006
        bb[ind] = 9528.5742
        ind = ((salt >= 34.131) & (salt < 34.29))
        mm[ind] = -127.2707
        bb[ind] = 4373.7895
        ind = ((salt >= 34.29) & (salt < 34.478))
        mm[ind] = 34.7556
        bb[ind] = -1182.0779
        ind = ((salt >= 34.478) & (salt < 35))
        mm[ind] = 401.7916
        bb[ind] = -13836.8132
        oxygen = mm*salt + bb
        # Limit values.
        oxygen[oxygen > 450] = 450
        oxygen[oxygen < 0] = 0
        return oxygen
    elif vn == ['TIC']:
        if np.nanmax(salt) > 36:
            print('Salt out of range for TIC regression')
        # Salinity vs. TIC [uM]
        # TIC = mm*salt + bb
        mm = 0*salt
        bb = 0*salt
        ind = (salt < 31.887)
        mm[ind] = 27.7967
        bb[ind] = 1112.2027
        ind = ((salt >= 31.887) & (salt < 33.926))
        mm[ind] = 147.002
        bb[ind] = -2688.8534
        ind = ((salt >= 33.926) & (salt < 34.197))
        mm[ind] = 352.9123
        bb[ind] = -9674.5448
        ind = ((salt >= 34.197) & (salt < 34.504))
        mm[ind] = 195.638
        bb[ind] = -4296.2223
        ind = ((salt >= 34.504) & (salt < 35))
        mm[ind] = -12.7457
        bb[ind] = 2893.77
        TIC = mm*salt + bb
        return TIC
    elif vn == ['alkalinity']:
        if np.nanmax(salt) > 36:
            print('Salt out of range for alkalinity regression')
        # Salinity vs. alkalinity [uM]
        # alkalinity = mm*salt + bb
        ind = (salt < 31.477)
        mm[ind] = 37.0543
        bb[ind] = 1031.0726
        ind = ((salt >= 31.477) & (salt < 33.915))
        mm[ind] = 48.5821
        bb[ind] = 668.2143
        ind = ((salt >= 33.915) & (salt < 35))
        mm[ind] = 246.2214
        bb[ind] = -6034.6841
        alkalinity = mm*salt + bb
        return alkalinity
        
def make_ini_file(nc_dir):
    # Initial condition, copied from first time of ocean_clm.nc
    ds1 = nc.Dataset(nc_dir + 'ocean_clm.nc', mode='r')
    # name output file
    ini_fn = nc_dir + 'ocean_ini.nc'
    # get rid of the old version, if it exists
    try:
        os.remove(ini_fn)
    except OSError:
        pass # assume error was because the file did not exist
    ds2 = nc.Dataset(ini_fn, 'w', format='NETCDF3_CLASSIC')
    # Copy dimensions
    for dname, the_dim in ds1.dimensions.items():
        if 'time' in dname:
            ds2.createDimension(dname, 1)
        else:
            ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
    # Copy variables
    for v_name, varin in ds1.variables.items():
        outVar = ds2.createVariable(v_name, varin.datatype, varin.dimensions)
        # Copy variable attributes, {} is a dict comprehension, cool!
        outVar.setncatts({k: varin.getncattr(k).replace('climatology','').strip() for k in varin.ncattrs()})
        if varin.ndim > 1:
            outVar[:] = varin[0,:]
        else:
            outVar[:] = varin[0]
    ds1.close()
    ds2.close()
    print('-Writing ocean_ini.nc')

def make_bry_file(nc_dir):
    # Boundary conditions, copied from edges of ocean_clm.nc
    ds1 = nc.Dataset(nc_dir + 'ocean_clm.nc', mode='r')
    # name output file
    bry_fn = nc_dir + 'ocean_bry.nc'
    # get rid of the old version, if it exists
    try:
        os.remove(bry_fn)
    except OSError:
        pass # assume error was because the file did not exist
    ds2 = nc.Dataset(bry_fn, 'w', format='NETCDF3_CLASSIC')
    # Copy dimensions
    for dname, the_dim in ds1.dimensions.items():
        ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
    # Copy parts of variables
    for v_name, varin in ds1.variables.items():
        if varin.ndim in [3,4]: 
            for bname in ['north', 'south', 'east', 'west']:
                outname = v_name + '_' + bname
                if bname in ['north', 'south']:
                    outdims = tuple([item for item in varin.dimensions if item[:3] != 'eta'])
                elif bname in ['west', 'east']:
                    outdims = tuple([item for item in varin.dimensions if item[:2] != 'xi'])
                outVar = ds2.createVariable(outname, varin.datatype, outdims)    
                outVar.setncatts({k: varin.getncattr(k).replace('climatology','').strip() for k in varin.ncattrs()})
                if varin.ndim == 4:
                    if bname == 'north':
                        outVar[:] = varin[:,:,-1,:]
                    elif bname == 'south':
                        outVar[:] = varin[:,:,0,:]
                    elif bname == 'east':
                        outVar[:] = varin[:,:,:,-1]
                    elif bname == 'west':
                        outVar[:] = varin[:,:,:,0]
                elif varin.ndim == 3:
                    if bname == 'north':
                        outVar[:] = varin[:,-1,:]
                    elif bname == 'south':
                        outVar[:] = varin[:,0,:]
                    elif bname == 'east':
                        outVar[:] = varin[:,:,-1]
                    elif bname == 'west':
                        outVar[:] = varin[:,:,0]
        else:
            outname = v_name
            outdims = tuple([item for item in varin.dimensions])
            outVar = ds2.createVariable(outname, varin.datatype, outdims)    
            outVar.setncatts({k: varin.getncattr(k).replace('climatology','').strip() for k in varin.ncattrs()})
            outVar[:] = varin[:]    
    ds1.close()
    ds2.close()
    print('-Writing ocean_bry.nc')
