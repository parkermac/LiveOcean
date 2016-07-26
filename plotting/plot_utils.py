# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 16:02:43 2016

@author: PM5

Module of basic utilities for plotting.  The goal is to make the code in
pfun less cumbersome to write and edit.

"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import zfun
import matfun

def auto_lims(fld):
    """
    A convenience function for automatically setting color limits.
    Input: a numpy array (masked OK)
    Output: tuple of good-guess colorsclale limits for a pcolormesh plot.
    """
    flo = max(fld.mean() - 3*fld.std(), np.nanmin(fld))
    fhi = min(fld.mean() + 3*fld.std(), np.nanmax(fld))
    return (flo, fhi)

def add_bathy_contours(ax, fn, depth_levs = []):
    if len(depth_levs) == 0:
        depth_levs = [100, 200, 500, 1000, 2000, 3000]
    ds = nc.Dataset(fn)
    h = ds['h'][:]
    lon = ds['lon_rho'][:]
    lat = ds['lat_rho'][:]
    ds.close()
    ax.contour(lon, lat, h, depth_levs, colors='g')

def add_map_field(ax, fn, varname, slev=-1, vlims=(), cmap='rainbow'):
    ds = nc.Dataset(fn)
    cmap = plt.get_cmap(name=cmap)
    if 'lon_rho' in ds[varname].coordinates:
        x = ds['lon_psi'][:]
        y = ds['lat_psi'][:]
        v = ds[varname][0, slev, 1:-1, 1:-1].squeeze()
    elif 'lon_v' in ds[varname].coordinates:
        x = ds['lon_u'][:]
        y = ds['lat_u'][:]
        v = ds[varname][0, slev, :, 1:-1].squeeze()
    elif 'lon_u' in ds[varname].coordinates:
        x = ds['lon_v'][:]
        y = ds['lat_v'][:]
        v = ds[varname][0, slev, 1:-1, :].squeeze()
    if len(vlims) == 0:
        vlims = auto_lims(v)
    ds.close()
    cs = ax.pcolormesh(x, y, v, vmin=vlims[0], vmax=vlims[1], cmap=cmap)
    return cs

def add_windstress_flower(ax, fn, t_scl=0.2, t_leglen =0.1):
    # ADD MEAN WINDSTRESS VECTOR
    # t_scl: scale windstress vector (smaller to get longer arrows)
    # t_leglen: # Pa for wind stress vector legend
    ds = nc.Dataset(fn)
    taux = ds['sustr'][:].squeeze()
    tauy = ds['svstr'][:].squeeze()
    tauxm = taux.mean()
    tauym = tauy.mean()
    ds.close()
    ax.quiver([.85, .85] , [.25, .25], [tauxm, tauxm], [tauym, tauym],
        units='y', scale=t_scl, scale_units='y', color='k',
        transform=ax.transAxes)
    tt = 1./np.sqrt(2)
    t_alpha = 0.3
    ax.quiver([.85, .85] , [.25, .25],
        t_leglen*np.array([0,tt,1,tt,0,-tt,-1,-tt]),
        t_leglen*np.array([1,tt,0,-tt,-1,-tt,0,tt]),
        units='y', scale=t_scl, scale_units='y', color='k', alpha=t_alpha,
        transform=ax.transAxes)
    ax.text(.85, .12,'Windstress',
        horizontalalignment='center', alpha=t_alpha, transform=ax.transAxes)
    ax.text(.85, .15, str(t_leglen) + ' Pa',
        horizontalalignment='center', alpha=t_alpha, transform=ax.transAxes)

def add_info(ax, fn):
    # put info on plot
    [T] = zfun.get_basic_info(fn, getG=False, getS=False)
    ax.text(.95, .07, T['tm'].strftime('%Y-%m-%d'),
        horizontalalignment='right', transform=ax.transAxes)
    ax.text(.95, .04, T['tm'].strftime('%H:%M') + ' UTC',
        horizontalalignment='right', transform=ax.transAxes)

def add_coast(ax, fn_coast=''):
    if len(fn_coast) == 0:
        fn_coast = ('/Users/PM5/Documents/LiveOcean_data' +
            '/coast/pnw_coast_combined.mat')
    cmat = matfun.loadmat(fn_coast)
    ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)

def get_aa(fn):
    ds = nc.Dataset(fn)
    x = ds['lon_psi'][0,:]
    y = ds['lat_psi'][:,0]
    ds.close()
    aa = [x[0], x[-1], y[0], y[-1]]
    return aa





