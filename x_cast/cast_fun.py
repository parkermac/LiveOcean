#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 16:25:26 2018

@author: pm7

module for cast extraction

"""

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import numpy as np
import zfun
import zrfun
import netCDF4 as nc

def get_cast(gridname, tag, ex_name, date_string, station, lon_str, lat_str):

    # get the dict Ldir
    Ldir = Lfun.Lstart(gridname, tag)
    Ldir['gtagex'] = Ldir['gtag'] + '_' + ex_name
    Ldir['date_string'] = date_string
    Ldir['station'] = station
    Ldir['lon_str'] = lon_str
    Ldir['lat_str'] = lat_str
    
    # make sure the output directory exists
    outdir0 = Ldir['LOo'] + 'cast/'
    Lfun.make_dir(outdir0)
    outdir = outdir0 + Ldir['gtagex'] + '/'
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
    fn = indir + 'ocean_his_0021.nc' # approx. noon local standard time
    
    if os.path.isfile(fn):
        pass
    else:
        print('Not found: ' + fn)
        return
    
    G = zrfun.get_basic_info(fn, only_G=True)
    S = zrfun.get_basic_info(fn, only_S=True)

    lon = G['lon_rho']
    lat = G['lat_rho']
    mask = G['mask_rho']
    xvec = lon[0,:].flatten()
    yvec = lat[:,0].flatten()
    
    i0, i1, frx = zfun.get_interpolant(np.array([float(Ldir['lon_str'])]), xvec)
    j0, j1, fry = zfun.get_interpolant(np.array([float(Ldir['lat_str'])]), yvec)
    i0 = int(i0)
    j0 = int(j0)
    # find indices of nearest good point
    if mask[j0,i0] == 1:
        print('- ' + station + ': point OK')
    elif mask[j0,i0] == 0:
        print('- ' + station + ':point masked')
        i0, j0 = get_ij_good(lon, lat, xvec, yvec, i0, j0, mask)
        new_lon = xvec[i0]
        new_lat = yvec[j0]
        Lon = np.array(new_lon)
        Lat = np.array(new_lat)
   
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
        
    print('-- Working on date: ' + dt)
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
        Ldir['station'] + '_' +
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
            
    foo.close()
    
def get_ij_good(lon, lat, xvec, yvec, i0, j0, mask):
    # find the nearest unmasked point
    #
    # starting point
    lon0 = lon[j0,i0]
    lat0 = lat[j0,i0]
    pad = 5 # how far to look (points)
    # indices of box to search over
    imax = len(xvec)-1
    jmax = len(yvec)-1
    I = np.arange(i0-pad, i0+pad)
    J = np.arange(j0-pad,j0+pad)
    # account for out-of-range points
    if I[0] < 0:
        I = I - I[0]
    if I[-1] > imax:
        I = I - (I[-1] - imax)
    if J[0] < 0:
        J = J - J[0]
    if J[-1] > jmax:
        J = J - (J[-1] - jmax)
    ii, jj = np.meshgrid(I, J)
    # sub arrays
    llon = lon[jj,ii]
    llat = lat[jj,ii]
    xxx, yyy = zfun.ll2xy(llon, llat, lon0, lat0)
    ddd = np.sqrt(xxx**2 + yyy**2) # distance from original point
    mmask = mask[jj,ii] 
    mm = mmask==1 # Boolean array of good points
    dddm = ddd[mm] # vector of good distances
    # indices of best point
    igood = ii[mm][dddm==dddm.min()][0]
    jgood = jj[mm][dddm==dddm.min()][0]
    #
    return igood, jgood

