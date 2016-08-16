# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 13:12:17 2016

@author: PM5
"""

def P_roms_sect(in_dict):
    # plots a section (distance, z)

    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE

    # GET DATA
    G, S, T = zfun.get_basic_info(fn)
    ds = nc.Dataset(fn,'r')
    h = G['h']
    zeta = ds.variables['zeta'][:].squeeze()
    zr = zfun.get_z(h, zeta, S, only_rho=True)

    varname = 'temp'

    sectvar = ds.variables[varname][:].squeeze()
    clims = {'salt':(28,35),'temp':(4,14)}

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

    # get coastline
    if len(fn_coast) != 0:
        # get the coastline
        cmat = matfun.loadmat(fn_coast)

    # CREATE THE SECTION
    # which track to use for the section
    if True:
        tracks_path = Ldir['data'] + 'tracks/'
        #which_track = 'lat47Track'
        which_track = 'jdf2psTrack'
        # get the track to interpolate onto
        mat = matfun.loadmat(tracks_path + which_track + '.mat')
        x = mat['x']
        y = mat['y']
    else:
        # create track by hand
        which_track = 'By hand'
        x = np.linspace(lon.min(), -124, 500)
        y = 48.5 * np.ones(x.shape)
    # create dist
    earth_rad = 6371e3 # m
    xrad = np.pi * x /180
    yrad = np.pi * y / 180
    dx = earth_rad * np.cos(yrad[1:]) * np.diff(xrad)
    dy = earth_rad * np.diff(yrad)
    ddist = np.sqrt(dx**2 + dy**2)
    dist = np.zeros(len(x))
    dist[1:] = ddist.cumsum()/1000 # km
    # find the index of zero
    it0 = zfun.get_interpolant(np.zeros(1), dist)
    idist0 = int(it0[0][0]) + int((it0[0][2].round()))
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
    xita = zfun.get_interpolant(x, xx)
    yita = zfun.get_interpolant(y, yy)
    # and prepare them to do the bilinear interpolation
    #xita = np.array(xit)
    #yita = np.array(yit)
    col0 = xita[:, 0].astype(int)
    col1 = xita[:, 1].astype(int)
    colf = xita[:, 2]
    colff = 1 - colf
    row0 = yita[:, 0].astype(int)
    row1 = yita[:, 1].astype(int)
    rowf = yita[:, 2]
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
            fldi = fldi.data # just the data, not the mask
        v3[fname] = fldi
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
    plt.close()
    fig = plt.figure(figsize=(20,10))

    # map
    ax = fig.add_subplot(131)
    cmap = plt.get_cmap(name='Set3')
    ax.contourf(lon, lat, zbot, 50, vmin=-3000, vmax=0, cmap=cmap)
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    ax.contour(lon, lat, zbot, [-2000, -1000, -500, -200, -100], colors='g')
    ax.set_xlim(lon.min(), lon.max())
    ax.set_ylim(lat.min(), lat.max())
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=10, markerfacecolor='w',
    markeredgecolor='r', markeredgewidth=2)
    # test the interpolation routine
    #ax.plot(v2['lon'], v2['lat'], '*k') # RESULT: it works
    if len(fn_coast) != 0:
        ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
    zfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Bathymetry and Section Track')

    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    #ax.set_ylim(v3['zr'].min(), 5)
    ax.set_ylim(-1000, 5)
    cmap = plt.get_cmap(name='rainbow')
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], v3['sectvarf'],  cmap=cmap)
    cs.set_clim(clims[varname])
    fig.colorbar(cs)
    cs = ax.contour(v3['distf'], v3['zrf'], v3['sectvarf'],
        np.arange(0, 40, .25), colors='k')
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title(which_track + ' ' + varname)

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
    else:
        plt.show()
    return out_dict

def P_tracks(in_dict):

    """
    Code to experiment with ways to create particle tracks.  Aimed at
    quick interactive visualization.
    """

    # START
    fig = plt.figure(figsize=figsize)
    ds = nc.Dataset(in_dict['fn'])
    vlims = in_dict['vlims'].copy()
    out_dict['vlims'] = vlims

    # PLOT CODE

    import netCDF4 as nc

    # GET DATA
    # run some code
    G, S, T = zfun.get_basic_info(fn)
    ds = nc.Dataset(fn,'r')
    u = ds.variables['u'][0, -1, :, :].squeeze()
    v = ds.variables['v'][0, -1, :, :].squeeze()
    #salt = ds.variables['salt'][0, -1, :, :].squeeze()
    chl = ds.variables['phytoplankton'][0, -1, :, :].squeeze()
    ds.close()
    # set masked values to 0
    ud = u.data; ud[G['mask_u']==False] = 0
    vd = v.data; vd[G['mask_v']==False] = 0

    # get coastline
    if len(fn_coast) != 0:
        # get the coastline
        cmat = matfun.loadmat(fn_coast)


    # make initial track locations
    aa = (G['lon_rho'][0,0],G['lon_rho'][0,-1],
        G['lat_rho'][0,0],G['lat_rho'][-1,0])
    # aa = [-124, -122.3, 47, 49] # Puget Sound Override
    daax = aa[1] - aa[0]
    daay = aa[3] - aa[2]
    mlat = np.mean(aa[2:])
    clat = np.cos(np.deg2rad(mlat))
    axrat = clat * daax / daay
    nngrid = 50
    nr = nngrid
    nc = round(nngrid * axrat)
    npts = nr * nc
    if False:
        #random grid
        x = np.random.uniform(aa[0], aa[1], npts)
        y = np.random.uniform(aa[2], aa[3], npts)
        #123.5 122
        #47 49

    else:
        xx = np.linspace(aa[0], aa[1], nc)
        yy = np.linspace(aa[2], aa[3], nr)
        XX,YY = np.meshgrid(xx,yy)
        x = XX.flatten()
        y = YY.flatten()

    # set track integration parameters
    dt = 1*3600. # time step (sec)
    nt = 24 # number of time steps
    RE = zfun.earth_rad(mlat)
    lonu = G['lon_u'][0, :]
    latu = G['lat_u'][:, 0]
    lonv = G['lon_v'][0, :]
    latv = G['lat_v'][:, 0]

    # initialize output arrays
    # the last time level of the locations will be nan, so
    # that we can plot tracks efficiently
    x2 = np.nan * np.ones((npts,nt+1))
    y2 = x2.copy()
    t2 = np.zeros((npts,nt+1))
    x2[:,0] = x
    y2[:,0] = y
    t2[:,0] = 0.

    for ii in range(1,nt):
        # get velocities at track locations
        ui = zfun.interp_scattered_on_plaid(x, y, lonu, latu, ud)
        vi = zfun.interp_scattered_on_plaid(x, y, lonv, latv, vd)
        # create distance step (m)
        dx = ui * dt
        dy = vi * dt
        # calculate new positions (degrees)
        dlon =  (dx/(clat*RE)) * (180./np.pi)
        dlat =  (dy/RE) * (180./np.pi)
        x = x + dlon
        y = y + dlat
        # store results
        x2[:,ii] = x
        y2[:,ii] = y
        t2[:,ii] = float(ii)

    # START PLOTTING
    plt.close()
    fig = plt.figure(figsize=(8,10))
    ax = fig.add_subplot(111)

    # PLOT FIELD
    cmap = plt.get_cmap(name='rainbow')
    v_lims = zfun.auto_lims(chl)
    cs = ax.pcolormesh(G['lon_psi'], G['lat_psi'], chl[1:-1,1:-1],
        vmin=v_lims[0], vmax=v_lims[1],  cmap = cmap)#, alpha=.7)
    ax.axis(aa)
    zfun.dar(ax)
    fig.colorbar(cs, ax=ax, extend='both')
    fs = 14
    ax.set_xlabel('Longitude', fontsize=fs)
    ax.set_ylabel('Latitude', fontsize=fs)

    # add coastline
    if len(fn_coast) != 0:
        ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)


    # PLOT TRACKS
    # plots lines with varying width
    from matplotlib.collections import LineCollection
    xx = x2.flatten()
    yy = y2.flatten()
    tt = t2.flatten()
    maxwidth = 2.
    lwidths = tt[:-1]*(maxwidth/nt) # drop last point
    points = np.array([xx, yy]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, linewidths=lwidths, colors='white')
    ax.add_collection(lc)

    # plot a dot
    #hr = np.mod(int(fn[-7:-3]),12)
    #ax.plot(x2[:,hr],y2[:,hr],'.r',markersize=2)

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
    else:
        plt.show()
    return out_dict