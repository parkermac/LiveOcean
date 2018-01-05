"""
Module of plotting functions.
"""

# imports
#
# because the calling function, pan_plot.py, has already put alpha on the
# path we assume it is on the path here too.
import numpy as np
import netCDF4 as nc

import Lfun
Ldir = Lfun.Lstart()
if Ldir['env'] == 'pm_mac': # mac version
    pass
elif Ldir['env'] == 'pm_fjord': # fjord version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt
import pickle

import cmocean as cmo

from importlib import reload
import zfun
import zrfun
import matfun
import pfun; reload(pfun)

# function for color limits
def get_in_dict(plot_type):
    #%% choices
    in_dict = dict()

    # COLOR LIMITS
    vlims = dict()
    # If you use () then the limits will be set by the first plot
    # and then held constant at those levels thereafter.    
    vlims['salt'] = (15, 35)#(28, 34)
    vlims['temp'] = (7, 18)
    vlims['NO3'] = (0, 40)
    vlims['phytoplankton'] = (0,30)#(0, 40)
    vlims['zooplankton'] = (0, 4)
    vlims['oxygen'] = (4, 8) # (0, 4) for bottom DO (ml L-1), (4, 8) for surface
    vlims['TIC'] = (2000, 2400) # (2000,2400) for surface
    vlims['alkalinity'] = (2000,2400)
    vlims['PH'] = (7, 8.5)#(6, 9)
    vlims['ARAG'] = (0, 3)#(0, 3)
    vlims['Ldetritus'] = ()
    # custom choices based on plot_type   
    if plot_type == 'P_layer':
        vlims['oxygen'] = (.5, 1) # (0, 4) for bottom DO (ml L-1), (4, 8) for surface
        vlims['TIC'] = (2450, 2650) # (2000,2400) for surface
    in_dict['vlims'] = vlims
        
    # OTHER
    in_dict['z_level'] = -250 # z level to plot
        
    return in_dict

# module defaults (available inside the methods)

# colormaps
cmap_dict = {'salt': 'rainbow', #cmo.cm.haline, 'gist_ncar'
             'temp': 'jet',#cmo.cm.thermal, #'bwr', cmo.cm.thermal, 'nipy_spectral'
             'NO3': cmo.cm.dense,
             'phytoplankton': 'jet',#cmo.cm.algae,
             'zooplankton': cmo.cm.matter,
             'oxygen': 'jet',#cmo.cm.oxy,
             'TIC': 'rainbow', #cmo.cm.matter,
             'alkalinity': 'rainbow', #cmo.cm.solar,
             'PH': 'jet',
             'ARAG': 'rainbow',
             'Ldetritus': 'rainbow'}

# units (after multiplying by fac)
units_dict = {'salt': '',
             'temp': ' $(^{\circ}C)$',
             'NO3': ' $(\mu mol\ L^{-1})$',
             'phytoplankton': ' $(\mu g\ chl\ C\ L^{-1})$',
             'zooplankton': ' $(\mu g\ chl\ C\ L^{-1})$',
             'oxygen': ' $(ml\ L^{-1})$',
             'TIC': ' $(\mu mol\ L^{-1})$',
             'alkalinity': ' $(\mu\ equivalents\ L^{-1})$',
             'PH': '',
             'ARAG': '',
             'Ldetritus': ''}

# scaling factors
fac_dict =  {'salt': 1,
             'temp': 1,
             'NO3': 1,
             'phytoplankton': 2.5,
             'zooplankton': 2.5,
             'oxygen': 0.032/1.42903, # convert mmol m-3 to ml L-1
             'TIC': 1,
             'alkalinity': 1,
             'PH': 1,
             'ARAG': 1,
             'Ldetritus': 1}
             
tstr_dict = {'salt': 'Salinity',
             'temp': 'Temperature',
             'NO3': 'Nitrate',
             'phytoplankton': 'Phytoplankton',
             'zooplankton': 'Zooplankton',
             'oxygen': 'DO',
             'TIC': 'DIC',
             'alkalinity': 'Alkalinity',
             'PH': 'pH',
             'ARAG': '$\Omega_{arag}$',
             'Ldetritus': 'Ldetritus'}

figsize = (13,8) # laptop
# figsize = (18,10) # big screen
out_dict = dict()

def P_basic(in_dict):
    # This creates, and optionally saves, a basic plot of surface fields
    # from a ROMS history file.
    # INPUT a dict containing:
    #   fn: text string with the full path name of the history file to plot
    #   fn_out: text string with full path of output file name
    #   in_dict: a tuple with optional information to pass to the plot
    # OUTPUT: either a screen image or a graphics file, and a dict of other
    # information such as axis limits.

    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE
    
    # HACKS
    auto_vlims = False
    #
    new_vlims = True
    if new_vlims==True:
        vlims['salt'] = (24,34)#(28, 34)
        vlims['temp'] = (5,11)
        
    # panel 1
    vn = 'salt'
    tstr = 'Surface ' + tstr_dict[vn]
    ax = fig.add_subplot(121)
    vn = 'salt'
    
    if auto_vlims:
        vlims[vn] = ()
    
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn])
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(tstr + units_dict[vn])
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)
    # panel 2
    ax = fig.add_subplot(122)
    vn = 'temp'
    
    if auto_vlims:
        vlims[vn] = ()
    
    tstr = 'Surface ' + tstr_dict[vn]
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn])
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(tstr + units_dict[vn])
    pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict

def P_salish(in_dict):
    # like basic, but the second panel focuses on the Salish Sea

    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE
    
    # HACKS
    auto_vlims = False
    #
    new_vlims = True
    if new_vlims==True:
        vlims['salt'] = (24,34)#(28, 34)
        vlims['temp'] = (5,11)
        
    # panel 1
    vn = 'salt'
    tstr = 'Surface ' + tstr_dict[vn]
    ax = fig.add_subplot(121)
    vn = 'salt'
    
    if auto_vlims:
        vlims[vn] = ()
    
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn])
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(tstr + units_dict[vn])
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)
    # panel 2
    ax = fig.add_subplot(122)
    vn = 'salt'
    
    if auto_vlims:
        vlims[vn] = ()
    
    tstr = 'Surface ' + tstr_dict[vn]
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=(26,32), cmap=cmap_dict[vn], fac=fac_dict[vn])
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis([-124, -122, 47, 49.5])
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(tstr + units_dict[vn])
    #pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict

def P_basic2D(in_dict):
    # For 2D fields.

    # START
    fig = plt.figure(figsize=(20,8))
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE
    # panel 1
    
    ax = fig.add_subplot(131)
    # we do this one by hand because zeta is not a masked array
    vn = 'zeta'
    x = ds['lon_psi'][:]
    y = ds['lat_psi'][:]
    v = ds[vn][0, 1:-1, 1:-1].squeeze()
    m = ds['mask_rho'][1:-1, 1:-1].squeeze()
    vm = np.ma.masked_where(m==0, v)
    vlims = pfun.auto_lims(vm)
    cs = ax.pcolormesh(x, y, vm, vmin=vlims[0], vmax=vlims[1], cmap='rainbow')
    out_dict['vlims'][vn] = vlims
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(vn)
    pfun.add_info(ax, in_dict['fn'])
    
    # panel 2
    ax = fig.add_subplot(132)
    vn = 'ubar'
    x = ds['lon_u'][:]
    y = ds['lat_u'][:]
    v = ds[vn][0, :, :].squeeze()
    m = ds['mask_u'][:, :].squeeze()
    vm = np.ma.masked_where(m==0, v)
    #vlims = pfun.auto_lims(vm)
    vlims = (vm.min(), vm.max())
    cs = ax.pcolormesh(x, y, vm, vmin=vlims[0], vmax=vlims[1], cmap='rainbow')
    out_dict['vlims'][vn] = vlims
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(vn)
    
    # panel 3
    ax = fig.add_subplot(133)
    vn = 'vbar'
    x = ds['lon_v'][:]
    y = ds['lat_v'][:]
    v = ds[vn][0, :, :].squeeze()
    m = ds['mask_v'][:, :].squeeze()
    vm = np.ma.masked_where(m==0, v)
    #vlims = pfun.auto_lims(vm)
    vlims = (vm.min(), vm.max())
    cs = ax.pcolormesh(x, y, vm, vmin=vlims[0], vmax=vlims[1], cmap='rainbow')
    out_dict['vlims'][vn] = vlims
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(vn)
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict

def P_basicN(in_dict):
    # Like P_basic, but optimized for the new nested grid sal0

    # START
    fig = plt.figure(figsize=(16,8))
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE
    #depth_levs = [50, 100, 150, 200]
    vlims['salt'] = (25,32)
    vlims['temp'] = (6,7.5)
    # panel 1
    vn = 'salt'
    tstr = 'Surface ' + tstr_dict[vn]
    ax = fig.add_subplot(121)
    vn = 'salt'
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn])
    fig.colorbar(cs)
    #pfun.add_bathy_contours(ax, ds, depth_levs=depth_levs)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(tstr + units_dict[vn])
    pfun.add_info(ax, in_dict['fn'])
    #pfun.add_windstress_flower(ax, ds, t_scl=1, t_leglen=0.1, center=(.25,.4))
    # panel 2
    ax = fig.add_subplot(122)
    vn = 'temp'
    tstr = 'Surface ' + tstr_dict[vn]
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap='rainbow', fac=fac_dict[vn])
    fig.colorbar(cs)
    #pfun.add_bathy_contours(ax, ds, depth_levs=depth_levs)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(tstr + units_dict[vn])
#    pfun.add_velocity_vectors(ax, ds, in_dict['fn'],
#                              v_scl=10, v_leglen=1, nngrid=70, center=(.25,.4))    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict

def P_dive_vort(in_dict):
    # plots surface fields of divergence and vorticity.
    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims
    # PLOT CODE
    u = ds['u'][0, -1, :, :]
    v = ds['v'][0, -1, :, :]
    G = zrfun.get_basic_info(in_dict['fn'], only_G=True)
    dive = ((np.diff(u, axis=1)/G['DX'][:, 1:-1])[1:-1, :]
            + (np.diff(v, axis = 0)/G['DY'][1:-1, :])[:, 1:-1])   
    x = G['lon_psi'] # matrix
    y = G['lat_psi'] # matrix
    dxp = zfun.interp2(x, y, G['lon_rho'], G['lat_rho'], G['DX'])
    dyp = zfun.interp2(x, y, G['lon_rho'], G['lat_rho'], G['DY'])
    vort = np.diff(v, axis=1)/dxp - np.diff(u, axis=0)/dyp
    
    aa = pfun.get_aa(ds)
    #aa = [-122.95, -122.55, 47.6, 48]
    scl = 1e-4
    # panel 1
    ax = fig.add_subplot(121)
    cs = plt.pcolormesh(G['lon_psi'], G['lat_psi'], dive, cmap='bwr',
                        vmin=-scl, vmax=scl)
    tstr = 'Surface Divergence (1/s)'
    #fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(tstr)
    pfun.add_info(ax, in_dict['fn'])
    
    # panel 2
    ax = fig.add_subplot(122)
    cs = plt.pcolormesh(G['lon_rho'], G['lat_rho'], vort, cmap='bwr',
                        vmin=-scl, vmax=scl)
    tstr = 'Surface Vorticity (1/s)'
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(tstr)

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict
    
def P_pH_Arag(in_dict):
    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE
    # panel 1
    ax = fig.add_subplot(121)
    vn = 'PH'
    tstr = 'Surface ' + tstr_dict[vn]
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn],
            do_mask_salish=True)
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(tstr + units_dict[vn])
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)
    # panel 2
    ax = fig.add_subplot(122)
    vn = 'ARAG'
    tstr = 'Surface ' + tstr_dict[vn]
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn],
            do_mask_salish=True)
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(tstr + units_dict[vn])
    pfun.add_velocity_vectors(ax, ds, in_dict['fn'])

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict

def P_Carbon(in_dict):
    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE
    # panel 1
    ax = fig.add_subplot(121)
    vn = 'TIC'
    tstr = 'Surface ' + tstr_dict[vn]
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn],
            do_mask_salish=True)
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(tstr + units_dict[vn])
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)
    # panel 2
    ax = fig.add_subplot(122)
    vn = 'alkalinity'
    tstr = 'Surface ' + tstr_dict[vn]
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn],
            do_mask_salish=True)
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(tstr + units_dict[vn])
    pfun.add_velocity_vectors(ax, ds, in_dict['fn'])

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict

def P_bio4(in_dict):
    # START
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims
    vn_list = ['NO3', 'phytoplankton', 'zooplankton', 'oxygen']
    NP = len(vn_list)
    if False:
        NR = np.maximum(1, np.ceil(np.sqrt(NP)).astype(int))
        NC = np.ceil(np.sqrt(NP)).astype(int)
        figsize = (16,16)
    else:
        NR = 1
        NC = NP
        figsize=(23, 7)
    fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=figsize,
                             squeeze=False)
    cc = 0
    for vn in vn_list:
        ir = int(np.floor(cc/NC))
        ic = int(cc - NC*ir)
        # PLOT CODE
        ax = axes[ir, ic]
        try:
            vlims[vn]
        except KeyError:
            vlims[vn] = ()
        if vn == 'oxygen':
            cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
                    vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn],
                    slev=0)
            ax.text(.95, .04, 'BOTTOM',
                    horizontalalignment='right', transform=ax.transAxes)
        else:
            cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
                    vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn])
        fig.colorbar(cs, ax=ax)
        if ic == 0:
            pfun.add_bathy_contours(ax, ds, txt=True)
        else:
            pfun.add_bathy_contours(ax, ds)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        if ir == NR-1:
            ax.set_xlabel('Longitude')
        if ic == 0:
            ax.set_ylabel('Latitude')
        ax.set_title(tstr_dict[vn] + units_dict[vn])
        if cc == 0:
            pfun.add_info(ax, in_dict['fn'])
        cc += 1

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict

def P_layer(in_dict):
    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])    
    vlims = in_dict['vlims'].copy()
    
    # set variables to plt
    #vn_list = ['Ldetritus','TIC']
    vn_list = ['salt', 'temp']
    #vn_list = ['NO3', 'temp']
    for vn in vn_list: # use auto scaling
        vlims[vn] = ()
    # and override
    vlims['Ldetritus'] = (0, 0.01)
    vlims['TIC'] = (2350, 2450)
    out_dict['vlims'] = vlims
    
    # PLOT CODE
    zfull = pfun.get_zfull(ds, in_dict['fn'], 'rho')
    
    # panel 1
    ax = fig.add_subplot(121)
    vn = vn_list[0]
    tstr = tstr_dict[vn]
    laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, in_dict['z_level'])
    cmap=cmap_dict[vn]
    fac=fac_dict[vn]
    if len(vlims[vn]) == 0:
        vlims[vn] = pfun.auto_lims(laym)
    out_dict['vlims'][vn] = vlims[vn]
    cs = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], fac*laym[1:-1,1:-1],
                       vmin=vlims[vn][0], vmax=vlims[vn][1], cmap=cmap)
    cb = fig.colorbar(cs)
    cb.formatter.set_useOffset(False)
    cb.update_ticks()
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(tstr + units_dict[vn] + ' on Z = ' + str(in_dict['z_level']) + ' m')
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)
    
    # panel 2
    ax = fig.add_subplot(122)
    vn = vn_list[1]
    tstr = tstr_dict[vn]
    laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, in_dict['z_level'])
    cmap=cmap_dict[vn]
    fac=fac_dict[vn]
    if len(vlims[vn]) == 0:
        vlims[vn] = pfun.auto_lims(laym)
    out_dict['vlims'][vn] = vlims[vn]
    cs = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], fac*laym[1:-1,1:-1],
                       vmin=vlims[vn][0], vmax=vlims[vn][1], cmap=cmap)
    cb = fig.colorbar(cs)
    cb.formatter.set_useOffset(False)
    cb.update_ticks()
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    # ****
    pfun.add_velocity_vectors(ax, ds, in_dict['fn'], v_scl=.5, v_leglen=0.1,
                              nngrid=120, zlev=in_dict['z_level'])
    # ****
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(tstr + units_dict[vn])

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict

def P_nest(in_dict):
    # plots a field nested inside the corresponding HYCOM field
    # - currently only works for salt and temp -

    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE
    # ** choose the variable to plot **
    vn = 'temp'
    # need to get Ldir, which means unpacking gtagex
    fn = in_dict['fn']
    gtagex = fn.split('/')[-3]
    gtx_list = gtagex.split('_')
    import Lfun
    Ldir = Lfun.Lstart(gtx_list[0], gtx_list[1])
    # Get HYCOM field
    hyvar_dict = {'temp':'t3d', 'salt':'s3d'}
    hvar = hyvar_dict[vn]
    f_string = fn.split('/')[-2]
    fnhy = (Ldir['LOo'] + Ldir['gtag'] + '/' + f_string +
        '/ocn/Data/' + hvar + '.nc')
    ds = nc.Dataset(fnhy)
    lon = ds.variables['lon'][:]
    lat = ds.variables['lat'][:]
    z = ds.variables['z'][:]
    # make sure zlev is in the HYCOM z
    zlev = zfun.find_nearest(z, in_dict['z_level'])
    nz = list(z).index(zlev)
    laym_hy = ds.variables[hvar + '_filt'][0,nz,:,:].squeeze()
    vlims = pfun.auto_lims(laym_hy)
    ds.close()
    aa = [lon.min(), lon.max(), lat.min(), lat.max()]
    # Get ROMS field
    ds = nc.Dataset(fn)
    zfull = pfun.get_zfull(ds, fn, 'rho')
    laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, zlev)

    # PLOTTING
    cmap = plt.get_cmap(name='gist_ncar')
    plt.close()
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=figsize, squeeze=False)
    # panel 1
    ax = axes[0,0]
    ax.set_title('HYCOM: z=' + str(zlev) + 'm, var=' + vn)
    cs = ax.pcolormesh(lon, lat, laym_hy,
                       vmin=vlims[0], vmax=vlims[1], cmap = cmap)
    ax.axis(aa)
    fig.colorbar(cs, ax=ax)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # panel 2
    ax = axes[0,1]
    ax.set_title('HYCOM + ROMS')
    cs = ax.pcolormesh(lon, lat, laym_hy,
                       vmin=vlims[0], vmax=vlims[1], cmap = cmap)
    ax.axis(aa)
    laymd = laym.data
    laymd[laym.mask] = 1e6
    ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], laymd[1:-1,1:-1],
                  vmin=vlims[0], vmax=vlims[1], cmap = cmap)
    fig.colorbar(cs, ax=ax)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict

def P_sect(in_dict):
    # plots a section (distance, z)

    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE

    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])

    vn = 'salt'
    try:
        vlims[vn]
    except KeyError:
        vlims[vn] = ()

    # CREATE THE SECTION
    # create track by hand
    if False:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -3500
        #x = np.linspace(lon.min(), -124, 500)
        if True:
            x = np.linspace(lon.min(), lon.max(), 500)
            y = 47 * np.ones(x.shape)
        else:
            y = np.linspace(lat.min(), lat.max(), 500)
            x = -126 * np.ones(y.shape)
    # or read one in
    else:
        import Lfun
        Ldir = Lfun.Lstart()
        tracks_path = Ldir['data'] + 'tracks_new/'
        which_track = 'HC_north.p'
        track_fn = tracks_path + which_track
        zdeep = -120
        # get the track to interpolate onto
        pdict = pickle.load(open(track_fn, 'rb'))
        xx = pdict['lon_poly']
        yy = pdict['lat_poly']
        for ii in range(len(xx)-1):
            x0 = xx[ii]
            x1 = xx[ii+1]
            y0 = yy[ii]
            y1 = yy[ii+1]
            nn = 20
            if ii == 0:
                x = np.linspace(x0, x1, nn)
                y = np.linspace(y0,y1, nn)
            else:
                x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
                y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))

    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING

    # panel 1
    ax = fig.add_subplot(1, 3, 1)
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=(24,32), cmap=cmap_dict[vn], fac=fac_dict[vn])
    #fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    #ax.axis(pfun.get_aa(ds))
    ax.axis([-123.25, -122.25, 47, 48.5])
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Bathymetry and Section Track')
    pfun.add_info(ax, in_dict['fn'])
    #pfun.add_windstress_flower(ax, ds)
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
    markeredgecolor='r', markeredgewidth=2)

    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    #vlims = pfun.auto_lims(v3['sectvarf'])
    vlims=(29,31)
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], v3['sectvarf'],
                       vmin=vlims[0], vmax=vlims[1], cmap=cmap_dict[vn])
    fig.colorbar(cs)
    cs = ax.contour(v3['distf'], v3['zrf'], v3['sectvarf'],
        np.linspace(np.floor(vlims[0]), np.ceil(vlims[1]), 20),
        colors='k', linewidths=0.5)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    tstr = tstr_dict[vn]
    ax.set_title(tstr)
    
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict
    
def P_sectA(in_dict):
    # plots a map and several sections
    # designed for analytical runs, like aestus1

    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims
    
    vn = 'salt'
    try:
        vlims[vn]
    except KeyError:
        vlims[vn] = ()
        
    # PLOTTING

    # map and section lines
    ax1 = fig.add_subplot(3,1,1)
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax1, ds, vn,
            vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn])
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax1, ds)
    pfun.add_coast(ax1)
    ax1.axis([-.5, 1, 44.8, 45.2])
    pfun.dar(ax1)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')
    pfun.add_info(ax1, in_dict['fn'], fs=9)
    
    # thalweg section
    x = np.linspace(-0.5, 1, 500)
    y = 45 * np.ones(x.shape)
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    ax = fig.add_subplot(3,1,2)
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(-25, 2)
    #vlims = pfun.auto_lims(v3['sectvarf'])
    vlims = vlims[vn]
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], v3['sectvarf'],
                       vmin=vlims[0], vmax=vlims[1], cmap=cmap_dict[vn])
    #fig.colorbar(cs)
    cs = ax.contour(v3['distf'], v3['zrf'], v3['sectvarf'],
        np.linspace(1, 35, 35), colors='k', linewidths=.5,)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title(vn)
    # add line to map plot
    ax1.plot(x, y, '-k', linewidth=2)
    ax1.plot(x[idist0], y[idist0], 'ok', markersize=10, markerfacecolor='w',
        markeredgecolor='k', markeredgewidth=2)

    for ii in range(3):
        # cross-sections
        y = np.linspace(44.9, 45.1, 200)
        x = 0.3*ii * np.ones(y.shape)
        v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
        # section
        ax = fig.add_subplot(3,3,ii+7)
        ax.plot(dist, v2['zbot'], '-k', linewidth=2)
        ax.plot(dist, v2['zeta'], '-b', linewidth=1)
        ax.set_xlim(0, 22.5)
        #ax.set_xlim(dist.min(), dist.max())
        ax.set_ylim(-25, 2)
        #vlims = pfun.auto_lims(v3['sectvarf'])
        cs = ax.pcolormesh(v3['distf'], v3['zrf'], v3['sectvarf'],
                           vmin=vlims[0], vmax=vlims[1], cmap=cmap_dict[vn])
        #fig.colorbar(cs)
        cs = ax.contour(v3['distf'], v3['zrf'], v3['sectvarf'],
            np.linspace(1, 35, 35), colors='k', linewidths=.5,)
        ax.set_xlabel('Distance (km)')
        if ii==0:
            ax.set_ylabel('Z (m)')
        # add line to map plot
        ax1.plot(x, y, '-k', linewidth=2)
        ax1.plot(x[idist0], y[idist0], 'ok', markersize=10, markerfacecolor='w',
            markeredgecolor='k', markeredgewidth=2)
            
    fig.tight_layout()
    

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict

def P_tracks_MERHAB(in_dict):
    # Use tracker to create surface drifter tracks for MERHAB 
    # It automatically makes tracks for as long as there are
    # hours in the folder, as a folder of movie frames.

    # START
    fig = plt.figure(figsize=(12, 8))
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # TRACKS
    import os
    import sys
    pth = os.path.abspath('../tracker')
    if pth not in sys.path:
        sys.path.append(pth)
    import trackfun
    import pickle

    # need to get Ldir, which means unpacking gtagex
    fn = in_dict['fn']
    gtagex = fn.split('/')[-3]
    gtx_list = gtagex.split('_')
    import Lfun
    Ldir = Lfun.Lstart(gtx_list[0], gtx_list[1])
    
    # make a list of history files in this folder,
    in_dir = fn[:fn.rindex('/')+1]
    fn_list_raw = os.listdir(in_dir)
    fn_list = []
    for item in fn_list_raw:
        if 'ocean_his' in item:
            fn_list.append(in_dir + item)
    fn_list.sort()
    
    # trim fn_list to end at the selected hour
    fn_list_full = fn_list.copy() # but save the full list
    fn_list = fn_list[:fn_list.index(fn)+1]
    # and estimate the number of days,
    ndays = round(len(fn_list)/24)
    # and use the CURRENT file for the map field overlay
    ds = nc.Dataset(in_dict['fn'])

    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    T0 = zrfun.get_basic_info(fn_list[0], only_T=True)

    if len(fn_list) == 2:
        # only do the tracking at the start
        
        # some run specifications
        ic_name = 'test' # 'jdf' or 'cr' or etc.
        dir_tag = 'forward' # 'forward' or 'reverse'
        method = 'rk4' # 'rk2' or 'rk4'
        surface = True # Boolean, True for trap to surface
        windage = 0.0 # a small number >= 0 [0.02]
        ndiv = 1 # number of divisions to make between saves for the integration

    
        # # Evenly spread over whole domain
        # x0 = G['lon_rho'][0, 1]
        # x1 = G['lon_rho'][0, -2]
        # y0 = G['lat_rho'][1, 0]
        # y1 = G['lat_rho'][-2, 0]
        # nyp = 30
        # mlr = np.pi*(np.mean([y0, y1]))/180
        # xyRatio = np.cos(mlr) * (x1 - x0) / (y1 - y0)
        # lonvec = np.linspace(x0, x1, (nyp * xyRatio).astype(int))
        # latvec = np.linspace(y0, y1, nyp)
        # lonmat, latmat = np.meshgrid(lonvec, latvec)
    
        if True:
            # standard MERHAB version
            nyp = 7
            x0 = -126
            x1 = -125
            y0 = 48
            y1 = 49
            mlr = np.pi*(np.mean([y0, y1]))/180
            xyRatio = np.cos(mlr) * (x1 - x0) / (y1 - y0)
            lonvec = np.linspace(x0, x1, (nyp * xyRatio).astype(int))
            latvec = np.linspace(y0, y1, nyp)
            lonmat_1, latmat_1 = np.meshgrid(lonvec, latvec)
            x0 = -125
            x1 = -124
            y0 = 44
            y1 = 45
            mlr = np.pi*(np.mean([y0, y1]))/180
            xyRatio = np.cos(mlr) * (x1 - x0) / (y1 - y0)
            lonvec = np.linspace(x0, x1, (nyp * xyRatio).astype(int))
            latvec = np.linspace(y0, y1, nyp)
            lonmat_2, latmat_2 = np.meshgrid(lonvec, latvec)
            lonmat = np.concatenate((lonmat_1.flatten(), lonmat_2.flatten()))
            latmat = np.concatenate((latmat_1.flatten(), latmat_2.flatten()))
        else:
            # new version with more release points
            nyp = 35
            x0 = -126
            x1 = -124
            y0 = 44
            y1 = 49
            mlr = np.pi*(np.mean([y0, y1]))/180
            xyRatio = np.cos(mlr) * (x1 - x0) / (y1 - y0)
            lonvec = np.linspace(x0, x1, (nyp * xyRatio).astype(int))
            latvec = np.linspace(y0, y1, nyp)
            lonmat, latmat = np.meshgrid(lonvec, latvec)

            
        plon00 = lonmat.flatten()
        plat00 = latmat.flatten()
        pcs00 = np.array([-.05]) # unimportant when surface=True

        # save some things in Ldir
        Ldir['gtagex'] = gtagex
        Ldir['ic_name'] = ic_name
        Ldir['dir_tag'] = dir_tag
        Ldir['method'] = method
        Ldir['surface'] = surface
        Ldir['windage'] = windage
        Ldir['ndiv'] = ndiv

        # make the full IC vectors, which will have equal length
        # (one value for each particle)
        NSP = len(pcs00)
        NXYP = len(plon00)
        plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        pcs0 = np.ones((NXYP,NSP)) * pcs00.reshape(1,NSP)
        plon0 = plon0.flatten()
        plat0 = plat0.flatten()
        pcs0 = pcs0.flatten()

        # DO THE TRACKING
        import time
        tt0 = time.time()
    
        P, Gtr, Str = trackfun.get_tracks(fn_list_full, plon0, plat0, pcs0,
                                      dir_tag, method, surface, ndiv, windage)
        print('  took %0.1f seconds' % (time.time() - tt0))
        
        fo = in_dict['fn_out']
        out_dir = fo[:fo.rindex('/')+1]
        out_fn = out_dir + 'tracks.p'
        pickle.dump(P, open(out_fn, 'wb'))
    else:
        fo = in_dict['fn_out']
        out_dir = fo[:fo.rindex('/')+1]
        out_fn = out_dir + 'tracks.p'
        P = pickle.load(open(out_fn, 'rb'))

    # PLOT CODE

    # panel 1
    ax = fig.add_subplot(121)
    vn = 'salt'
    tstr = 'Surface ' + tstr_dict[vn] +' and ' + str(ndays) + ' day Tracks'
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn], alpha = .5,
            do_mask_salish=False)
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    fs1 = 16
    ax.set_xlabel('Longitude', fontsize=fs1)
    ax.set_ylabel('Latitude', fontsize=fs1)
    ax.set_title(tstr, fontsize=fs1)
    
    fs = fs1 - 6
    ax.text(.98, .10, T0['tm'].strftime('%Y-%m-%d %H:%M'),
        horizontalalignment='right' , verticalalignment='bottom',
        transform=ax.transAxes, fontsize=fs)
    ax.text(.98, .07, 'to ' + T['tm'].strftime('%Y-%m-%d %H:%M'),
        horizontalalignment='right', verticalalignment='bottom',
        transform=ax.transAxes, fontsize=fs)
    ax.text(.98, .04, 'UTC',
        horizontalalignment='right', verticalalignment='bottom',
        transform=ax.transAxes, fontsize=fs)
    ax.text(.06, .04, fn.split('/')[-3],
        verticalalignment='bottom', transform=ax.transAxes,
        rotation='vertical', fontsize=fs)

    # add the tracks
    c_start = 'w'; s_start = 4
    c_end = 'r'; s_end = 6
    ntt = len(fn_list) -1
    ax.plot(P['lon'][:ntt,:], P['lat'][:ntt,:], '-k', linewidth=2)
    ax.plot(P['lon'][0,:],P['lat'][0,:],'o'+c_start,
            markersize=s_start, alpha = 1, markeredgecolor='k')
    ax.plot(P['lon'][ntt,:],P['lat'][ntt,:],'o'+c_end,
            markersize=s_end, alpha = 1, markeredgecolor='k')
    
    # add info about the tracks    
    x0 = .8; x1 = .9
    y0 = .25; y1 = .3
    ax.plot([x0, x1], [y0, y1], '-k', linewidth=2, transform=ax.transAxes)
    ax.plot(x0, y0,'o'+c_start,
            markersize=s_start, alpha = 1, markeredgecolor='k', transform=ax.transAxes)
    ax.plot(x1, y1,'o'+c_end,
            markersize=s_end, alpha = 1, markeredgecolor='k', transform=ax.transAxes)
    # and some labels
    ax.text(x0+.02, y0, 'start', horizontalalignment='left',
            verticalalignment='center', fontstyle='italic', transform=ax.transAxes)
    ax.text(x1-.02, y1, 'end', horizontalalignment='right',
            verticalalignment='center', fontstyle='italic', transform=ax.transAxes)
   
    ax = fig.add_subplot(122)
    vn = 'phytoplankton'
    tstr = 'Surface ' + tstr_dict[vn]
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
           vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn],
           do_mask_salish=True)
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=False)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    fs1 = 16
    ax.set_xlabel('Longitude', fontsize=fs1)
    #ax.set_ylabel('Latitude', fontsize=fs1)
    ax.set_title(tstr, fontsize=fs1)
   
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict
    
def P_aestus(in_dict):
    # designed for the analytical estuary-shelf runs

    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE
    # panel 1
    ax = fig.add_subplot(111)
    vn = 'salt'
    tstr = 'Surface ' + tstr_dict[vn]
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=(0,35), cmap=cmap_dict[vn], fac=fac_dict[vn])
    fig.colorbar(cs)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(tstr + units_dict[vn])
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_velocity_vectors(ax, ds, in_dict['fn'], nngrid=40)

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict

