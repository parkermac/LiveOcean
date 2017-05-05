# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 16:02:43 2016

@author: PM5

Module of basic utilities for plotting.  The goal is to make the code in
pfun less cumbersome to write and edit.

"""

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
import zrfun

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import pandas as pd

def topfig():
    """
    Positions the figure at the top left of the screen, and puts it
    on top of the Spyder window.  Copied from Stack Overflow.
    """
    try:
        figmgr = plt.get_current_fig_manager()
        figmgr.canvas.manager.window.raise_()
        geom = figmgr.window.geometry()
        x,y,dx,dy = geom.getRect()
        figmgr.window.setGeometry(10, 10, dx, dy)
    except AttributeError:
        # this only works, and is only needed, in Spyder
        # so working from the terminal we just skip it.
        pass

def dar(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.
    """
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.sin(np.pi*yav/180))

def add_coast(ax, dir0=Ldir['data']):
    fn = dir0 + 'coast/coast_pnw.p'
    C = pd.read_pickle(fn)
    ax.plot(C['lon'].values, C['lat'].values, '-k', linewidth=0.5)

def get_coast(dir0=Ldir['data']):
    fn = dir0 + 'coast/coast_pnw.p'
    C = pd.read_pickle(fn)
    lon = C['lon'].values
    lat = C['lat'].values
    return lon, lat

def auto_lims(fld):
    """
    A convenience function for automatically setting color limits.
    Input: a numpy array (masked OK)
    Output: tuple of good-guess colorsclale limits for a pcolormesh plot.
    """
    flo = np.nanmax([np.nanmean(fld) - 3*np.nanstd(fld), np.nanmin(fld)])
    fhi = np.nanmin([np.nanmean(fld) + 3*np.nanstd(fld), np.nanmax(fld)])
    return (flo, fhi)

def get_units(ds, vn):
    try:
        units = ds[vn].units
    except AttributeError:
        units = ''
    return units

def add_bathy_contours(ax, ds, depth_levs = [], txt=False):
    # this should work with ds being a history file Dataset, or the G dict.
    h = ds['h'][:]
    lon = ds['lon_rho'][:]
    lat = ds['lat_rho'][:]
    c1 = 'k'
    c2 = 'k'
    if len(depth_levs) == 0:
        ax.contour(lon, lat, h, [200], colors=c1, linewidths=0.5)
        ax.contour(lon, lat, h, [2000], colors=c2, linewidths=0.5)
        if txt==True:
            ax.text(.95, .95, '200 m', color=c1,
                    horizontalalignment='right',
                    transform=ax.transAxes)
            ax.text(.95, .92, '2000 m', color=c2,
                    horizontalalignment='right',transform=ax.transAxes)
    else:
        ax.contour(lon, lat, h, depth_levs, colors='g')

def add_map_field(ax, ds, varname, slev=-1, vlims=(), cmap='rainbow',
                  fac=1, alpha=1, do_mask_salish=False):
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
    
    if do_mask_salish:
        v = mask_salish(v, ds['lon_rho'][1:-1, 1:-1], ds['lat_rho'][1:-1, 1:-1])  
    
    cs = ax.pcolormesh(x, y, v*fac, vmin=vlims[0], vmax=vlims[1], cmap=cmap, alpha=alpha)
    return cs, vlims

def add_velocity_vectors(ax, ds, fn, v_scl=3, v_leglen=0.5, nngrid=80):
    # v_scl: scale velocity vector (smaller to get longer arrows)
    # v_leglen: m/s for velocity vector legend
    # GET DATA
    G = zrfun.get_basic_info(fn, only_G=True)
    u = ds['u'][0, -1, :, :].squeeze()
    v = ds['v'][0, -1, :, :].squeeze()
    # ADD VELOCITY VECTORS
    # set masked values to 0
    ud = u.data; ud[G['mask_u']==False] = 0
    vd = v.data; vd[G['mask_v']==False] = 0
    # create interpolant
    import scipy.interpolate as intp
    ui = intp.interp2d(G['lon_u'][0, :], G['lat_u'][:, 0], ud)
    vi = intp.interp2d(G['lon_v'][0, :], G['lat_v'][:, 0], vd)
    # create regular grid
    aaa = ax.axis()
    daax = aaa[1] - aaa[0]
    daay = aaa[3] - aaa[2]
    axrat = np.cos(np.deg2rad(aaa[2])) * daax / daay
    x = np.linspace(aaa[0], aaa[1], round(nngrid * axrat))
    y = np.linspace(aaa[2], aaa[3], nngrid)
    xx, yy = np.meshgrid(x, y)
    # interpolate to regular grid
    uu = ui(x, y)
    vv = vi(x, y)
    mask = uu != 0
    # plot velocity vectors
    ax.quiver(xx[mask], yy[mask], uu[mask], vv[mask],
        units='y', scale=v_scl, scale_units='y', color='k')
    ax.quiver([.7, .7] , [.05, .05], [v_leglen, v_leglen],
              [v_leglen, v_leglen],
        units='y', scale=v_scl, scale_units='y', color='k',
        transform=ax.transAxes)
    ax.text(.75, .05, str(v_leglen) + ' $ms^{-1}$',
        horizontalalignment='left', transform=ax.transAxes)

def add_windstress_flower(ax, ds, t_scl=0.2, t_leglen=0.1):
    # ADD MEAN WINDSTRESS VECTOR
    # t_scl: scale windstress vector (smaller to get longer arrows)
    # t_leglen: # Pa for wind stress vector legend
    taux = ds['sustr'][:].squeeze()
    tauy = ds['svstr'][:].squeeze()
    tauxm = taux.mean()
    tauym = tauy.mean()
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

def add_info(ax, fn, fs=12):
    # put info on plot
    T = zrfun.get_basic_info(fn, only_T=True)
    ax.text(.95, .07, T['tm'].strftime('%Y-%m-%d'),
        horizontalalignment='right', transform=ax.transAxes, fontsize=fs)
    ax.text(.95, .04, T['tm'].strftime('%H:%M') + ' UTC',
        horizontalalignment='right', transform=ax.transAxes, fontsize=fs)
    ax.text(.06, .04, fn.split('/')[-3],
        verticalalignment='bottom', transform=ax.transAxes,
        rotation='vertical', fontsize=fs)

def get_aa(ds):
    x = ds['lon_psi'][0,:]
    y = ds['lat_psi'][:,0]
    aa = [x[0], x[-1], y[0], y[-1]]
    return aa
    
def get_aa_ex(ds):
    x = ds['lon_psi_ex'][0,:]
    y = ds['lat_psi_ex'][:,0]
    aa = [x[0], x[-1], y[0], y[-1]]
    return aa

def get_zfull(ds, fn, which_grid):
    # get zfull field on "which_grid" ('rho', 'u', or 'v')
    G, S, T = zrfun.get_basic_info(fn)
    zeta = 0 * ds.variables['zeta'][:].squeeze()
    zr_mid = zrfun.get_z(G['h'], zeta, S, only_rho=True)
    zr_bot = -G['h'].reshape(1, G['M'], G['L']).copy()
    zr_top = zeta.reshape(1, G['M'], G['L']).copy()
    zfull0 = make_full((zr_bot, zr_mid, zr_top))
    if which_grid == 'rho':
        zfull = zfull0
    elif which_grid == 'u':
        zfull = zfull0[:, :, 0:-1] + np.diff(zfull0, axis=2)/2
    elif which_grid == 'v':
        zfull = zfull0[:, 0:-1, :] + np.diff(zfull0, axis=1)/2
    return zfull

def get_laym(ds, zfull, mask, vn, zlev):
    # make the layer
    fld_mid = ds[vn][:].squeeze()
    fld = make_full((fld_mid,))
    zlev_a = zlev * np.ones(1)
    lay = get_layer(fld, zfull, zlev_a)
    lay[mask == False] = np.nan
    laym = np.ma.masked_where(np.isnan(lay), lay)
    return laym

def get_layer(fld, zfull, which_z):
    """
    Creates a horizontal slice through a 3D ROMS data field.  It is very fast
    because of the use of "choose"
    Input:
        fld (3D ndarray) of the data field to slice
        z (3D ndarray) of z values (like from make_full)
        which_z (ndarray of length 1) of the z value for the layer
    Output:
        lay (2D ndarray) fld on z == which_z,
            with np.nan where it is not defined
    """
    N, M, L = fld.shape # updates N for full fields
    Nmax = 30
    ii = np.arange(0,N,Nmax)
    ii = np.concatenate((ii,np.array([N,])))
    fld0 = np.nan * np.zeros((M, L), dtype=int)
    fld1 = np.nan * np.zeros((M, L), dtype=int)
    z0 = np.nan * np.zeros((M, L), dtype=int)
    z1 = np.nan * np.zeros((M, L), dtype=int)
    # NOTE: need fewer than 32 layers to use "choose"
    # so we split the operation into steps in this loop
    j = 0
    while j < len(ii)-1:
        i_lo = ii[j]
        i_hi = min(ii[j+1] + 1, ii[-1]) # overlap by 1
        NN = i_hi - i_lo # the number of levels in this chunk
        this_zr = zfull[i_lo:i_hi].copy()
        this_fld = fld[i_lo:i_hi].copy()
        zm = this_zr < which_z
        ind0 = np.zeros((M, L), dtype=int)
        ind1 = np.zeros((M, L), dtype=int)
        ind0 = (zm == True).sum(0) - 1 # index of points below which_z
        ind1 = ind0 + 1 # index of points above which_z
        # dealing with out-of-bounds issues
        # note 0 <= ind0 <= NN-1
        # and  1 <= ind1 <= NN
        # make ind1 = ind0 for out of bounds cases
        ind0[ind0 == -1] = 0 # fix bottom case
        ind1[ind1 == NN] = NN-1 # fix top case
        # and now cells that should be masked have equal indices
        this_mask = ind0 != ind1
        this_fld0 = ind0.choose(this_fld)
        this_fld1 = ind1.choose(this_fld)
        this_z0 = ind0.choose(this_zr)
        this_z1 = ind1.choose(this_zr)
        fld0[this_mask] = this_fld0[this_mask]
        fld1[this_mask] = this_fld1[this_mask]
        z0[this_mask] = this_z0[this_mask]
        z1[this_mask] = this_z1[this_mask]
        j += 1
    # do the interpolation
    dz = z1 - z0
    dzf = which_z - z0
    dz[dz == 0] = np.nan
    fr = dzf / dz
    lay = fld0*(1 - fr) + fld1*fr
    return lay

def make_full(flt):
    """
    Adds top and bottom layers to array fld. This is intended for 3D ROMS data
    fields that are on the vertical rho grid, and where we want (typically for
    plotting purposes) to extend this in a smart way to the sea floor and the
    sea surface.
    NOTE: input is always a tuple.  If just sending a single array pack it
    as zfun.make_full((arr,))
    Input:
        flt is a tuple with either 1 ndarray (fld_mid,),
        or 3 ndarrays (fld_bot, fld_mid, fld_top)
    Output: fld is the "full" field
    """
    if len(flt)==3:
       fld = np.concatenate(flt, axis=0)
    elif len(flt)==1:
        if len(flt[0].shape) == 3:
            fld_mid = flt[0]
            N, M, L = fld_mid.shape
            fld_bot = fld_mid[0].copy()
            fld_bot = fld_bot.reshape(1, M, L).copy()
            fld_top = fld_mid[-1].copy()
            fld_top = fld_top.reshape(1, M, L).copy()
            fld = np.concatenate((fld_bot, fld_mid, fld_top), axis=0)
        elif len(flt[0].shape) == 2:
            fld_mid = flt[0]
            N, M = fld_mid.shape
            fld_bot = fld_mid[0].copy()
            fld_bot = fld_bot.reshape(1, M).copy()
            fld_top = fld_mid[-1].copy()
            fld_top = fld_top.reshape(1, M).copy()
            fld = np.concatenate((fld_bot, fld_mid, fld_top), axis=0)
        elif len(flt[0].shape) == 1:
            fld_mid = flt[0]
            fld_bot = fld_mid[0].copy()
            fld_top = fld_mid[-1].copy()
            fld = np.concatenate((fld_bot, fld_mid, fld_top), axis=0)
    return fld
    
def mask_salish(fld, lon, lat):
    """
    Mask out map fields inside the Salish Sea.   
    Input:
        2D fields of data (masked array), and associated lon and lat
        all must be the same shap
    Output:
        The data field, now masked in the Salish Sea.
    """
    x = [-125.5, -123.5, -122, -122]
    y = [50, 46.8, 46.8, 50]
    V = np.ones((len(x),2))
    V[:,0] = x
    V[:,1] = y
    P = mpath.Path(V)
    Rlon = lon.flatten()
    Rlat = lat.flatten()
    R = np.ones((len(Rlon),2))
    R[:,0] = Rlon
    R[:,1] = Rlat
    RR = P.contains_points(R) # boolean    
    fld = np.ma.masked_where(RR.reshape(lon.shape), fld)
    return fld





