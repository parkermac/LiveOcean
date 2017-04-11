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
    vlims['salt'] = (28, 34)
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
    in_dict['z_level'] = -1500 # z level to plot
        
    return in_dict

# module defaults (available inside the methods)

# colormaps
cmap_dict = {'salt': 'jet', #cmo.cm.haline,
             'temp': 'jet',#'bwr', #cmo.cm.thermal,
             'NO3': cmo.cm.dense,
             'phytoplankton': 'jet',#cmo.cm.algae,
             'zooplankton': cmo.cm.matter,
             'oxygen': 'jet',#cmo.cm.oxy,
             'TIC': cmo.cm.matter,
             'alkalinity': cmo.cm.solar,
             'PH': 'jet',
             'ARAG': 'rainbow',
             'Ldetritus': 'jet'}

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

figsize = (12,8) # laptop
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
    # panel 1
    vn = 'salt'
    tstr = 'Surface ' + tstr_dict[vn]
    ax = fig.add_subplot(121)
    vn = 'salt'
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
    [G] = zrfun.get_basic_info(in_dict['fn'], getG=True, getS=False, getT=False)
    dive = (np.diff(u, axis=1)/G['DX'][:, 1:-1])[1:-1, :] + (np.diff(v, axis = 0)/G['DY'][1:-1, :])[:, 1:-1]
    
    x = G['lon_psi'] # matrix
    y = G['lat_psi'] # matrix
    vv = zfun.interp_scattered_on_plaid(x, y, G['lon_rho'][0,:], G['lat_rho'][:,0], G['DX'])
    dxp = np.reshape(vv, G['lon_psi'].shape)
    vv = zfun.interp_scattered_on_plaid(x, y, G['lon_rho'][0,:], G['lat_rho'][:,0], G['DY'])
    dyp = np.reshape(vv, G['lat_psi'].shape)
    vort = np.diff(v, axis=1)/dxp - np.diff(u, axis=0)/dyp
    
    aa = [-122.95, -122.55, 47.6, 48]
    scl = 5e-4
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

def P_bio2(in_dict):
    # START
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims
    vn_list = ['phytoplankton', 'oxygen']
    NP = len(vn_list)
    NR = 1
    NC = NP
    #figsize = (18,10)
    fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=figsize,
                             squeeze=False)
    cc = 0
    for vn in vn_list:
        ir = int(np.floor(cc/NC))
        ic = int(cc - NC*ir)
        # PLOT CODE
        ax = axes[ir, ic]
        tstr = tstr_dict[vn]
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
        ax.set_title(tstr + units_dict[vn])
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
    vlims['salt'] = ()
    vlims['temp'] = ()
    
    out_dict['vlims'] = vlims
    
    # PLOT CODE
    zfull = pfun.get_zfull(ds, in_dict['fn'], 'rho')
    
    # panel 1
    ax = fig.add_subplot(121)
    vn = 'salt'
    tstr = tstr_dict[vn]
    laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, in_dict['z_level'])
    cmap=cmap_dict[vn]
    fac=fac_dict[vn]
    #vlims[vn] = ()
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
    vn = 'temp'
    tstr = tstr_dict[vn]
    laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, in_dict['z_level'])
    cmap=cmap_dict[vn]
    fac=fac_dict[vn]
    #vlims[vn] = ()
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
    fig = plt.figure(figsize=(20,8))
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE
    from warnings import filterwarnings
    filterwarnings('ignore') # skip this warning message:
    #/Applications/anaconda/lib/python3.5/site-packages/
    #matplotlib/colors.py:581: RuntimeWarning:
    #invalid value encountered in less
    #cbook._putmask(xa, xa < 0.0, -1)

    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    h = G['h']
    zeta = ds['zeta'][:].squeeze()
    zr = zrfun.get_z(h, zeta, S, only_rho=True)

    vn = 'salt'
    try:
        vlims[vn]
    except KeyError:
        vlims[vn] = ()

    sectvar = ds[vn][:].squeeze()

    L = G['L']
    M = G['M']
    N = S['N']

    lon = G['lon_rho']
    lat = G['lat_rho']
    mask = G['mask_rho']
    maskr = mask.reshape(1, M, L).copy()
    mask3 = np.tile(maskr, [N, 1, 1])
    zbot = -h # don't need .copy() because of the minus operation

    # make sure fields are masked
    zeta[mask==False] = np.nan
    zbot[mask==False] = np.nan
    sectvar[mask3==False] = np.nan

    # CREATE THE SECTION
    # create track by hand
    if False:
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
        tracks_path = Ldir['data'] + 'tracks/'
        which_track = 'jdf2psTrack'
        zdeep = -400
        # get the track to interpolate onto
        mat = matfun.loadmat(tracks_path + which_track + '.mat')
        x = mat['x']
        y = mat['y']

    # create dist
    earth_rad = zfun.earth_rad(np.mean(lat[:,0])) # m
    xrad = np.pi * x /180
    yrad = np.pi * y / 180
    dx = earth_rad * np.cos(yrad[1:]) * np.diff(xrad)
    dy = earth_rad * np.diff(yrad)
    ddist = np.sqrt(dx**2 + dy**2)
    dist = np.zeros(len(x))
    dist[1:] = ddist.cumsum()/1000 # km
    # find the index of zero
    i0, i1, fr = zfun.get_interpolant(np.zeros(1), dist)
    idist0 = i0
    distr = dist.reshape(1, len(dist)).copy()
    dista = np.tile(distr, [N, 1]) # array
    # pack fields to process in dicts
    d2 = dict()
    d2['zbot'] = zbot
    d2['zeta'] = zeta
    d2['lon'] = lon
    d2['lat'] = lat
    d3 = dict()
    d3['zr'] = zr
    d3['sectvar'] = sectvar
    # get vectors describing the (plaid) grid
    xx = lon[1,:]
    yy = lat[:,1]
    col0, col1, colf = zfun.get_interpolant(x, xx)
    row0, row1, rowf = zfun.get_interpolant(y, yy)
    # and prepare them to do the bilinear interpolation
    colff = 1 - colf
    rowff = 1 - rowf
    # now actually do the interpolation
    # 2-D fields
    v2 = dict()
    for fname in d2.keys():
        fld = d2[fname]
        fldi = (rowff*(colff*fld[row0, col0] + colf*fld[row0, col1])
        + rowf*(colff*fld[row1, col0] + colf*fld[row1, col1]))
        if type(fldi) == np.ma.core.MaskedArray:
            fldi = fldi.data # just the data, not the mask
        v2[fname] = fldi
    # 3-D fields
    v3 = dict()
    for fname in d3.keys():
        fld = d3[fname]
        fldi = (rowff*(colff*fld[:, row0, col0] + colf*fld[:, row0, col1])
        + rowf*(colff*fld[:, row1, col0] + colf*fld[:, row1, col1]))
        if type(fldi) == np.ma.core.MaskedArray:
            fldid = fldi.data # just the data, not the mask
            fldid[fldi.mask == True] = np.nan
        v3[fname] = fldid
    v3['dist'] = dista # distance in km
    # make "full" fields by padding top and bottom
    nana = np.nan * np.ones((N + 2, len(dist))) # blank array
    v3['zrf'] = nana.copy()
    v3['zrf'][0,:] = v2['zbot']
    v3['zrf'][1:-1,:] = v3['zr']
    v3['zrf'][-1,:] = v2['zeta']
    #
    v3['sectvarf'] = nana.copy()
    v3['sectvarf'][0,:] = v3['sectvar'][0,:]
    v3['sectvarf'][1:-1,:] = v3['sectvar']
    v3['sectvarf'][-1,:] = v3['sectvar'][-1,:]
    #
    v3['distf'] = nana.copy()
    v3['distf'][0,:] = v3['dist'][0,:]
    v3['distf'][1:-1,:] = v3['dist']
    v3['distf'][-1,:] = v3['dist'][-1,:]
    # NOTE: should make this a function
    # (but note that it is specific to each variable)

    # PLOTTING

    # panel 1
    ax = fig.add_subplot(1, 3, 1)
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn])
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Bathymetry and Section Track')
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=10, markerfacecolor='w',
    markeredgecolor='r', markeredgewidth=2)

    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    vlims = pfun.auto_lims(v3['sectvarf'])
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], v3['sectvarf'],
                       vmin=vlims[0],
                       vmax=vlims[1],
                       cmap='rainbow')
    fig.colorbar(cs)
    cs = ax.contour(v3['distf'], v3['zrf'], v3['sectvarf'],
        np.linspace(np.floor(vlims[0]), np.ceil(vlims[1]), 20), colors='k')
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    tstr = tstr_dict[vn]
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
    
def P_sectA(in_dict):
    # plots a section (distance, z)
    # designed for analytical runs, like aestus1

    # START
    fig = plt.figure(figsize=(14,10))
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE
    from warnings import filterwarnings
    filterwarnings('ignore') # skip a warning message

    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    h = G['h']
    zeta = ds['zeta'][:].squeeze()
    zr = zrfun.get_z(h, zeta, S, only_rho=True)

    vn = 'salt'
    try:
        vlims[vn]
    except KeyError:
        vlims[vn] = ()

    sectvar = ds[vn][:].squeeze()

    L = G['L']
    M = G['M']
    N = S['N']

    lon = G['lon_rho']
    lat = G['lat_rho']
    mask = G['mask_rho']
    maskr = mask.reshape(1, M, L).copy()
    mask3 = np.tile(maskr, [N, 1, 1])
    zbot = -h # don't need .copy() because of the minus operation

    # make sure fields are masked
    zeta[mask==False] = np.nan
    zbot[mask==False] = np.nan
    sectvar[mask3==False] = np.nan

    # CREATE THE SECTION
    # create track by hand
    x = np.linspace(-0.5, 1, 500)
    y = 45 * np.ones(x.shape)

    # create dist
    earth_rad = zfun.earth_rad(np.mean(lat[:,0])) # m
    xrad = np.pi * x /180
    yrad = np.pi * y / 180
    dx = earth_rad * np.cos(yrad[1:]) * np.diff(xrad)
    dy = earth_rad * np.diff(yrad)
    ddist = np.sqrt(dx**2 + dy**2)
    dist = np.zeros(len(x))
    dist[1:] = ddist.cumsum()/1000 # km
    # find the index of zero
    i0, i1, fr = zfun.get_interpolant(np.zeros(1), dist)
    idist0 = i0
    distr = dist.reshape(1, len(dist)).copy()
    dista = np.tile(distr, [N, 1]) # array
    # pack fields to process in dicts
    d2 = dict()
    d2['zbot'] = zbot
    d2['zeta'] = zeta
    d2['lon'] = lon
    d2['lat'] = lat
    d3 = dict()
    d3['zr'] = zr
    d3['sectvar'] = sectvar
    # get vectors describing the (plaid) grid
    xx = lon[1,:]
    yy = lat[:,1]
    col0, col1, colf = zfun.get_interpolant(x, xx)
    row0, row1, rowf = zfun.get_interpolant(y, yy)
    # and prepare them to do the bilinear interpolation
    colff = 1 - colf
    rowff = 1 - rowf
    # now actually do the interpolation
    # 2-D fields
    v2 = dict()
    for fname in d2.keys():
        fld = d2[fname]
        fldi = (rowff*(colff*fld[row0, col0] + colf*fld[row0, col1])
        + rowf*(colff*fld[row1, col0] + colf*fld[row1, col1]))
        if type(fldi) == np.ma.core.MaskedArray:
            fldi = fldi.data # just the data, not the mask
        v2[fname] = fldi
    # 3-D fields
    v3 = dict()
    for fname in d3.keys():
        fld = d3[fname]
        fldi = (rowff*(colff*fld[:, row0, col0] + colf*fld[:, row0, col1])
        + rowf*(colff*fld[:, row1, col0] + colf*fld[:, row1, col1]))
        if type(fldi) == np.ma.core.MaskedArray:
            fldid = fldi.data # just the data, not the mask
            fldid[fldi.mask == True] = np.nan
        v3[fname] = fldid
    v3['dist'] = dista # distance in km
    # make "full" fields by padding top and bottom
    nana = np.nan * np.ones((N + 2, len(dist))) # blank array
    v3['zrf'] = nana.copy()
    v3['zrf'][0,:] = v2['zbot']
    v3['zrf'][1:-1,:] = v3['zr']
    v3['zrf'][-1,:] = v2['zeta']
    #
    v3['sectvarf'] = nana.copy()
    v3['sectvarf'][0,:] = v3['sectvar'][0,:]
    v3['sectvarf'][1:-1,:] = v3['sectvar']
    v3['sectvarf'][-1,:] = v3['sectvar'][-1,:]
    #
    v3['distf'] = nana.copy()
    v3['distf'][0,:] = v3['dist'][0,:]
    v3['distf'][1:-1,:] = v3['dist']
    v3['distf'][-1,:] = v3['dist'][-1,:]
    # NOTE: should make this a function
    # (but note that it is specific to each variable)

    # PLOTTING

    # panel 1
    ax = fig.add_subplot(211)
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=(10,35.5), cmap='rainbow', fac=fac_dict[vn])
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis([-.5, 1, 44.8, 45.2])
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Section Track')
    pfun.add_info(ax, in_dict['fn'])
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=10, markerfacecolor='w',
    markeredgecolor='r', markeredgewidth=2)

    # section
    ax = fig.add_subplot(212)
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(-25, 2)
    vlims = pfun.auto_lims(v3['sectvarf'])
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], v3['sectvarf'],
                       vmin=10,
                       vmax=35.5,
                       cmap='rainbow')
    fig.colorbar(cs)
    cs = ax.contour(v3['distf'], v3['zrf'], v3['sectvarf'],
        np.linspace(1, 34, 34), colors='k')
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title(varname)

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        pfun.topfig()
    return out_dict

def P_tracks(in_dict):
    # use tracker to create surface drifter tracks
    # to run this in ipython do something like
    # run pan_plot.py -d 2016.08.25
    # (a day that has a 72 hour forecast in it)
    # and then select the number for List type = test
    # and then select the number for Plot type = P_tracks
    # and it automatically makes tracks for as long as there are
    # hours in the folder

    # START
    fig = plt.figure(figsize=(12, 12))
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # TRACKS
    import os
    import sys
    pth = os.path.abspath('../tracker')
    if pth not in sys.path:
        sys.path.append(pth)
    import trackfun

    # need to get Ldir, which means unpacking gtagex
    fn = in_dict['fn']
    gtagex = fn.split('/')[-3]
    gtx_list = gtagex.split('_')
    import Lfun
    Ldir = Lfun.Lstart(gtx_list[0], gtx_list[1])

    # some run specifications
    ic_name = 'test' # 'jdf' or 'cr' or etc.
    dir_tag = 'forward' # 'forward' or 'reverse'
    method = 'rk4' # 'rk2' or 'rk4'
    surface = True # Boolean, True for trap to surface
    windage = 0.0 # a small number >= 0 [0.02]
    ndiv = 1 # number of divisions to make between saves for the integration

    G, S, T = zrfun.get_basic_info(fn)

    in_dir = fn[:fn.rindex('/')+1]
    fn_list_raw = os.listdir(in_dir)
    fn_list = []
    for item in fn_list_raw:
        if 'ocean_his' in item:
            fn_list.append(in_dir + item)
    ndays = round(len(fn_list)/24)
    
    x0 = G['lon_rho'][0, 1]
    x1 = G['lon_rho'][0, -2]
    y0 = G['lat_rho'][1, 0]
    y1 = G['lat_rho'][-2, 0]
    nyp = 30
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
    # run some code
    
    P, Gtr, Str = trackfun.get_tracks(fn_list, plon0, plat0, pcs0,
                                  dir_tag, method, surface, ndiv, windage)
    print('  took %0.1f seconds' % (time.time() - tt0))

    # PLOT CODE

    # panel 1
    ax = fig.add_subplot(111)
    vn = 'phytoplankton'
    tstr = 'Surface ' + tstr_dict[vn] +' and ' + str(ndays) + ' day Tracks'
    cs, out_dict['vlims'][vn] = pfun.add_map_field(ax, ds, vn,
            vlims=vlims[vn], cmap=cmap_dict[vn], fac=fac_dict[vn], alpha = .3)
    fig.colorbar(cs)
    #pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(tstr)
    pfun.add_info(ax, in_dict['fn'])

    # add the tracks
    ax.plot(P['lon'], P['lat'], '-k')
    ax.plot(P['lon'][0,:],P['lat'][0,:],'ok',
            markersize=6, alpha = 1, markeredgecolor='k')

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
    fig = plt.figure(figsize=(14,8))
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

