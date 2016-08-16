"""
Module of plotting functions.
"""

# imports
#
# because the calling function, pan_plot.py, has already put alpha on the
# path we assume it is on the path here too.
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

from importlib import reload
import zfun; reload(zfun)
import pfun; reload(pfun)

# module defaults (available inside the methods)

figsize = (18,10)
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
    # panel 1
    t_str = 'Surface Salinity'
    ax = fig.add_subplot(121)
    vn = 'salt'
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap='rainbow')
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(t_str)
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)
    # panel 2
    t_str = 'Surface Temperature'
    ax = fig.add_subplot(122)
    vn = 'temp'
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap='bwr')
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(t_str + ' (' + pfun.get_units(ds, vn) + ')')
    pfun.add_velocity_vectors(ax, ds, in_dict['fn'])

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    return out_dict

def P_carbon(in_dict):
    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE
    # panel 1
    t_str = 'Surface TIC'
    ax = fig.add_subplot(121)
    vn = 'TIC'
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap='YlOrBr')
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(t_str + ' (' + pfun.get_units(ds, vn) + ')')
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)
    # panel 2
    t_str = 'Surface Alkalinity'
    ax = fig.add_subplot(122)
    vn = 'alkalinity'
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap='RdYlGn')
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(t_str + ' (' + pfun.get_units(ds, vn) + ')')
    pfun.add_velocity_vectors(ax, ds, in_dict['fn'])

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    return out_dict

def P_bio(in_dict):
    # START
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    vn_list = ['NO3', 'phytoplankton', 'zooplankton', 'oxygen']
    cmap_list = ['Oranges', 'Greens', 'BuPu', 'Spectral']
    cmap_dict = dict(zip(vn_list, cmap_list))
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
        t_str = vn
        try:
            vlims[vn]
        except KeyError:
            vlims[vn] = ()
        if vn == 'oxygen':
            cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
                vlims=vlims[vn], cmap=cmap_dict[vn], slev=0)
            ax.text(.95, .04, 'BOTTOM',
                    horizontalalignment='right', transform=ax.transAxes)
        else:
            cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
                vlims=vlims[vn], cmap=cmap_dict[vn])
        fig.colorbar(cs, ax=ax)
        pfun.add_bathy_contours(ax, ds)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title(t_str + ' (' + pfun.get_units(ds, vn) + ')')
        if ir == NR-1:
            ax.set_xlabel('Longitude')
        if ic == 0:
            ax.set_ylabel('Latitude')
        ax.set_title(t_str + ' (' + pfun.get_units(ds, vn) + ')')
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
    return out_dict

def P_layer(in_dict):
    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE
    zfull = pfun.get_zfull(ds, in_dict['fn'], 'rho')
    # panel 1
    t_str = 'Salinity'
    ax = fig.add_subplot(121)
    vn = 'salt'
    laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, in_dict['z_level'])
    if len(vlims[vn]) == 0:
        vlims[vn] = pfun.auto_lims(laym)
    out_dict['vlims'][vn] = vlims[vn]
    cs = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], laym[1:-1,1:-1],
                       vmin=vlims[vn][0], vmax=vlims[vn][1], cmap='rainbow')
    cb = fig.colorbar(cs)
    cb.formatter.set_useOffset(False)
    cb.update_ticks()
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(t_str + ' on Z = ' + str(in_dict['z_level']) + ' m')
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)
    # panel 2
    t_str = 'Temperature'
    ax = fig.add_subplot(122)
    vn = 'temp'
    laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, in_dict['z_level'])
    if len(vlims[vn]) == 0:
        vlims[vn] = pfun.auto_lims(laym)
    out_dict['vlims'][vn] = vlims[vn]
    cs = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], laym[1:-1,1:-1],
                       vmin=vlims[vn][0], vmax=vlims[vn][1], cmap='bwr')
    cb = fig.colorbar(cs)
    cb.formatter.set_useOffset(False)
    cb.update_ticks()
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(t_str + ' (' + pfun.get_units(ds, vn) + ')')

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    return out_dict

def P_nest_plot(in_dict):
    # plots a field nested inside the corresponding HYCOM field
    # - currently only works for salt and temp -

    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE
    # ** choose the variable to plot **
    which_var = 'temp'
    # need to get Ldir, which means unpacking gtagex
    fn = in_dict['fn']
    gtagex = fn.split('/')[-3]
    gtx_list = gtagex.split('_')
    import Lfun
    Ldir = Lfun.Lstart(gtx_list[0], gtx_list[1])
    # Get HYCOM field
    hyvar_dict = {'temp':'t3d', 'salt':'s3d'}
    hvar = hyvar_dict[which_var]
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
    laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], which_var, zlev)

    # PLOTTING
    cmap = plt.get_cmap(name='gist_ncar')
    plt.close()
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=figsize, squeeze=False)
    # panel 1
    ax = axes[0,0]
    ax.set_title('HYCOM: z=' + str(zlev) + 'm, var=' + which_var)
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
    return out_dict
