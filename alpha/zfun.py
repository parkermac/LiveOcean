"""
This module contains functions for extracting data from ROMS history files.

Parker MacCready
"""
def interp_scattered_on_plaid(x, y, xvec, yvec, u):
    """
    Gets values of the field u at locations (x,y).

    All inputs and outputs are numpy arrays.

    Field u is defined on a plaid grid defined by vectors xvec and yvec.

    Locations (x,y) are defined by vectors whose elements specify
    individual points.  They are not a plaid grid, but can be
    scattered arbitrarily throughout the domain.

    Returns a vector ui the same length at x and y.

    Note that because it relies on "get_interpolant_fast" out-of-bounds
    values default to the edge values of field u.
    """
    # get interpolants
    xint = get_interpolant_fast(x,xvec)
    yint = get_interpolant_fast(y,yvec)

    # and use these to get the interpolated u and v
    xi0 = xint[:,0].astype(int)
    yi0 = yint[:,0].astype(int)
    xi1 = xint[:,1].astype(int)
    yi1 = yint[:,1].astype(int)
    xf = xint[:,2]
    yf = yint[:,2]

    # bi linear interpolation
    u00 = u[yi0,xi0]
    u10 = u[yi1,xi0]
    u01 = u[yi0,xi1]
    u11 = u[yi1,xi1]
    ui = (1-yf)*((1-xf)*u00 + xf*u01) + yf*((1-xf)*u10 + xf*u11)

    return ui

def get_interpolant_fast(x, xvec):
    """
    Returns info to allow fast interpolation.

    Input: data point(s) x and coordinate vector xvec
    (both must be 1-D numpy arrays)

    Output: indices into xvec that surround x,
    and the fraction 'a' into that segment to find x
    which I call "interpolants"
    returned as a 3-column numpy array with columns:
    [index below, index above, fraction]

    If the input is ON a point in xvec the default is to return
    the index of that point and the one above.
    """
    import numpy as np

    # input error checking
    if isinstance(x, np.ndarray) and isinstance(xvec, np.ndarray):
        pass # input type is good
    else:
        print('WARNING from get_interpolant_fast(): ' +
              'Inputs must be numpy arrays')

    # some preconditioning of the input
    x = x.flatten()
    xvec = xvec.flatten()

    # more error checking
    if not np.all(np.diff(xvec) > 0):
        print('WARNING from get_interpolant_fast(): ' +
            'xvec must be monotonic and increasing')

    nx = len(x)
    nxvec = len(xvec)

    X = x.reshape(nx, 1) # column vector
    xvec = xvec.reshape(1, nxvec)
    XVEC = xvec.repeat(nx, axis=0) # matrix

    # preallocate results array
    itp = np.zeros((nx,3))

    # calculate index columns
    mask = X >= XVEC
    # the above line broadcasts correctly even if nx = nxvec
    # because we forced X to be a column vector
    itp[:,0] = mask.sum(axis=1) - 1
    # these masks are used to handle values of x beyond the limits
    # of xvec
    lomask = itp[:,0] < 0
    himask = itp[:,0] > nxvec - 2
    itp[lomask, 0] = 0
    itp[himask, 0] = nxvec - 2
    itp[:,1] = itp[:,0] + 1

    xvec0 = xvec[0,itp[:,0].astype(int)]
    xvec1 = xvec[0,itp[:,1].astype(int)]
    frac = (x - xvec0)/(xvec1 - xvec0)
    itp[:,2] = frac
    itp[lomask, 2] = 0
    itp[himask, 2] = 1

    return itp


def get_interpolant(x, xvec):
    """
    OBSOLETE: should be phased out (in moor_0.py and pfun.py).
    
    Returns info to allow fast interpolation (I hope).

    Input: data point(s) x and coordinate vector xvec
    (both must be 1-D numpy arrays)

    Output: indices into xvec that surround x,
    and the fraction 'a' into that segment to find x
    which I call "interpolants"
    returned as a list of three-element tuples,
    with each tuple containing (index below, index above, fraction)
    """
    import numpy as np

    # input error checking (could also use "isinstance")
    if type(x) != np.ndarray or type(xvec) != np.ndarray:
        print('WARNING from get_interpolant(): Inputs must be numpy arrays')
        ind0 = ind1 = a = np.nan
        return zip(ind0, ind1, a)
    if not np.all(np.diff(xvec) > 0):
        print('WARNING from get_interpolant():' +
              'xvec must be monotonic and increasing')
        ind0 = ind1 = a = np.nan
        return zip(ind0, ind1, a)

    # some preconditioning of the input
    x = x.flatten()
    xvec = xvec.flatten()
    # xvec.sort()  # not needed because we check above

    # make array of indices
    ind = np.arange(len(xvec))

    # preallocate results vectors
    N = len(x)
    ind0 = np.zeros(N, dtype=int)
    ind1 = np.zeros(N, dtype=int)
    a = np.zeros(N, dtype=float)

    # calculate results
    n = 0 # a counter
    for xx in x: # surely we could find a faster way to to this!
        # calculate indices, with some choices about edge and out-of-bounds values
        if xx <= xvec[0]:
            ind0[n] = 0; ind1[n] = 1; a[n] = 0.
        elif xx >= xvec[-1]:
            ind0[n] = len(xvec) - 1; ind1[n] = len(xvec); a[n] = 1.
        else:
            mask = xvec < xx
            ind0[n] = ind[mask][-1]
            ind1[n] = ind0[n] + 1
            # calculate fraction
            dx = xvec[ind1[n]] - xvec[ind0[n]]
            dxp = xx - xvec[ind0[n]]
            a[n] = dxp / dx

        n += 1

    return zip(ind0, ind1, a)

def get_basic_info(fn, getG=True, getS=True, getT=True):
    """
    Gets grid, vertical coordinate, and time info from a ROMS NetCDF
    history file with full name 'fn'

    Input: the filename (with path if needed)

    Output: dicts G, S, and T

    Example calls:

    for more than one the dout list is unpacked automatically
    G, S, T = zfun.get_basic_info(fn)

    for getting just one use [] to just get the dict
    [T] = zfun.get_basic_info(fn, getG=False, getS=False)

    """

    import netCDF4 as nc
    import numpy as np

    ds = nc.Dataset(fn,'r')

    dout = [] # initialize an output list

    if getG:
        # get grid and bathymetry info
        g_varlist = ['h', 'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v',
        'lon_psi', 'lat_psi', 'mask_rho', 'mask_u', 'mask_v', 'pm', 'pn',]
        G = dict()
        for vv in g_varlist:
            G[vv] = ds.variables[vv][:]
        G['DX'] = 1/G['pm']
        G['DY'] = 1/G['pn']
        G['M'], G['L'] = np.shape(G['lon_rho']) # M = rows, L = columns
        # make the masks boolean
        G['mask_rho'] = G['mask_rho'] == 1
        G['mask_u'] = G['mask_u'] == 1
        G['mask_v'] = G['mask_v'] == 1
        dout.append(G)

    if getS:
        # get vertical sigma-coordinate info (vectors are bottom to top)
        s_varlist = ['s_rho', 's_w', 'hc', 'Cs_r', 'Cs_w', 'Vtransform']
        S = dict()
        for vv in s_varlist:
            S[vv] = ds.variables[vv][:]
        S['N'] = len(S['s_rho']) # number of vertical levels
        dout.append(S)

    if getT:
        # get time info
        t_varlist = ['ocean_time', 'dstart']
        T = dict()
        for vv in t_varlist:
            T[vv] = ds.variables[vv][:]
        #
        # find  time reference
        dstart = ds.variables['dstart']
        tu = dstart.units
        import re
        isdash = [m.start() for m in re.finditer('-', tu)]
        iscolon = [m.start() for m in re.finditer(':', tu)]
        #
        year = int(tu[isdash[0]-4:isdash[0]])
        month = int(tu[isdash[1]-2:isdash[1]])
        day = int(tu[isdash[1]+1:isdash[1]+3])
        #
        hour = int(tu[iscolon[0]-2:iscolon[0]])
        minute = int(tu[iscolon[1]-2:iscolon[1]])
        second = int(tu[iscolon[1]+1:iscolon[1]+3])
        #
        import datetime
        tt = datetime.datetime(year, month, day, hour, minute, second)
        delta = datetime.timedelta(0, int(T['ocean_time']))
        T['tm0'] = tt
        T['tm'] = tt + delta
        dout.append(T)

    return dout

def get_z(h, zeta, S, only_rho=False, only_w=False):
    """
    Used to calculate the z position of fields in a ROMS history file

    Input: arrays h (bathymetry depth) and zeta (sea surface height)
    which must be the same size, and dict S created by get_basic_info()

    Output: 3-D arrays of z_rho and z_w

    NOTE: one foible is that if you input arrays of h and zeta that are
    vectors of length VL, the output array (e.g. z_rho) will have size (N, VL)
    (i.e. it will never return an array with size (N, VL, 1), even if (VL, 1) was
    the input shape).  This is a result of the initial and final squeeze calls.
    """

    import numpy as np

    # input error checking (seems clumsy)
    if type(h) != np.ndarray or type(zeta) not in [np.ndarray, np.ma.core.MaskedArray] or type(S) != dict:
        print('WARNING from get_z(): Inputs must be numpy arrays')

    # number of vertical levels
    N = S['N']

    # remove singleton dimensions
    h = h.squeeze()
    zeta = zeta.squeeze()
    # ensure that we have enough dimensions
    h = np.atleast_2d(h)
    zeta = np.atleast_2d(zeta)
    # check that the dimensions are the same
    if h.shape != zeta.shape:
        print('WARNING from get_z(): h and zeta must be the same shape')
    M, L = h.shape

    if not only_w:
        # rho
        # create some useful arrays
        csr = S['Cs_r']
        csrr = csr.reshape(N, 1, 1).copy()
        Cs_r = np.tile(csrr, [1, M, L])
        H_r = np.tile(h.reshape(1, M, L).copy(), [N, 1, 1])
        Zeta_r = np.tile(zeta.reshape(1, M, L).copy(), [N, 1, 1])
        if S['hc'] == 0: # if hc = 0 the transform is simpler (and faster)
            z_rho = H_r*Cs_r + Zeta_r + Zeta_r*Cs_r
        elif S['hc'] != 0: # need to calculate a few more useful arrays
            sr = S['s_rho']
            srr = sr.reshape(N, 1, 1).copy()
            S_rho = np.tile(srr, [1, M, L])
            Hc_r = np.tile(S['hc'], [N, M, L])
            if S['Vtransform'] == 1:
                zr0 = (S_rho - Cs_r) * Hc_r + Cs_r*H_r
                z_rho = zr0 + Zeta_r * (1 + zr0/H_r)
            elif S['Vtransform'] == 2:
                zr0 = (S_rho*Hc_r + Cs_r*H_r) / (Hc_r + H_r)
                z_rho = Zeta_r + (Zeta_r + H_r)*zr0
        z_rho = z_rho.squeeze()

    if not only_rho:
        # w
        # create some useful arrays
        csw = S['Cs_w']
        csww = csw.reshape(N+1, 1, 1).copy()
        Cs_w = np.tile(csww, [1, M, L])
        H_w = np.tile(h.reshape(1, M, L).copy(), [N+1, 1, 1])
        Zeta_w = np.tile(zeta.reshape(1, M, L).copy(), [N+1, 1, 1])
        if S['hc'] == 0: # if hc = 0 the transform is simpler (and faster)
            z_w = H_w*Cs_w + Zeta_w + Zeta_w*Cs_w
        elif S['hc'] != 0: # need to calculate a few more useful arrays
            sw = S['s_w']
            sww = sw.reshape(N+1, 1, 1).copy()
            S_w = np.tile(sww, [1, M, L])    #
            Hc_w = np.tile(S['hc'], [N+1, M, L])
            if S['Vtransform'] == 1:
                zw0 = (S_w - Cs_w) * Hc_w + Cs_w*H_w
                z_w = zw0 + Zeta_w * (1 + zw0/H_w)
            elif S['Vtransform'] == 2:
                zw0 = (S_w*Hc_w  + Cs_w*H_w) / (Hc_w + H_w)
                z_w = Zeta_w + (Zeta_w + H_w)*zw0
        z_w = z_w.squeeze()

    # return results
    if (not only_rho) and (not only_w):
        return z_rho, z_w
    elif only_rho and (not only_w):
        return z_rho
    elif (not only_rho) and only_w:
        return z_w

def dar(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.

    Input: axes object

    Output: none (but it alters the plot)
    """
    import numpy as np
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.sin(np.pi*yav/180))

def labxy(ax, pos='ll'):
    """
    returns x and y for axis labels
    """
    aa = ax.axis()
    daax = aa[1] - aa[0]
    daay = aa[3] - aa[2]
    delx = 0.05
    dely = 0.05
    if pos == 'll':
        xt = aa[0] + delx * daax
        yt = aa[2] + dely * daay
    elif pos == 'lr':
        xt = aa[1] - delx * daax
        yt = aa[2] + dely * daay
    elif pos == 'ul':
        xt = aa[0] + delx * daax
        yt = aa[3] - dely * daay
    elif pos == 'ur':
        xt = aa[1] - delx * daax
        yt = aa[3] - dely * daay
    else:
        pass

    return xt, yt


def get_layer(fld, zr, which_z):
    """
    Creates a horizontal slice through a 3D ROMS data field.  It is very fast
    because of the use of "choose"

    Input:
        fld (3D ndarray) of the data field to slice
        z (3D ndarray) of z values
        which_z (ndarray of length 1) of the z value for the layer

    Output:
        lay (2D ndarray) fld on z == which_z, with np.nan where it is not defined
    """
    import numpy as np
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

        this_zr = zr[i_lo:i_hi].copy()
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

    Output:
        fld is the "full" field
    """
    import numpy as np

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

def interpolate4D(ds0, ds1, varname, itp, ics, ilat, ilon):
    """
    Linear interpolation in t, z, y, x, for finding a value
    inside a single gridbox in ROMS history files.

    Could rewrite to be more general as two functions, the first to
    get a 4D array, and the second to interpolate within it.

    Input:
        two NetCDF Datasets, at different times
        the name of the variable to work on
        four interpolant tuples, created by get_interpolant

    Output:
        the interpolated value of the variable
        should return nan if any elements are nan

    Need to test carefully,maybe using an analytic function.
    """
    import numpy as np

    box0 = ds0.variables[varname][0,ics[:2],ilat[:2],ilon[:2]].squeeze()
    box1 = ds1.variables[varname][0,ics[:2],ilat[:2],ilon[:2]].squeeze()

    at = np.array([itp[:2]])
    az = np.array([1-ics[2], ics[2]])
    ay = np.array([1-ilat[2], ilat[2]]).reshape((1,2))
    ax = np.array([1-ilon[2], ilon[2]]).reshape((1,1,2))

    bb0 =  (az*( ( ay*((ax*box0).sum(-1)) ).sum(-1) )).sum()
    bb1 =  (az*( ( ay*((ax*box1).sum(-1)) ).sum(-1) )).sum()
    ival = (at*np.array([bb0, bb1])).sum()

    return ival

def filt_AB8d(data):
    """
    % 8/28/2013  Parker MacCready
    % ** use ONLY with hourly data! **
    %
    % This applies the Austin-Barth 8 day filter to a single vector.
    %
    % NOTE this may be different than their definition - it just returns a
    % weighted average of the data over the previous 8 days from time t,
    % with the weighting decaying from 1 to 1/e at t - 8 days.  There are 192
    % hours in 8 days.

    Input:
        assumed to be a 1D numpy array

    Output:
        a vector of the same size you started with,
        padded with NaN's at the ends.
    """
    import numpy as np

    fl = 8*24;

    filt = np.exp(np.linspace(-1,0,fl))
    filt = filt/filt.sum();

    smooth = np.nan*data

    for ii in range(fl+1, len(data)):

        smooth[ii] = (filt*data[ii-fl+1: ii + 1]).sum()

    return smooth

def filt_hanning(data, n=40):
    """
    2005 (approx.) Jonathan Lilly, modified by Parker MacCready
    BUG: seems to require  n to be even
    Input:
        data, assumed to be a 1D numpy array
        n is the filter length, default = 40
        if n=1 it just returns matrix you gave it

    Output:
        a vector of the same size you started with,

    % by default it uses a Hanning window
    """
    import numpy as np

    if n == 1:
        smooth = data
    else:
        # make a Hanning window from scratch
        ff = np.cos(np.linspace(-np.pi,np.pi,n+2))[1:-1]
        filt = (1 + ff)/2
        filt = filt / filt.sum()
        a = round(n/2.)
        smooth = np.nan * data
        smooth[a: -a + 1]= np.convolve(data, filt, mode = 'valid')
        #smooth[:n+1] = np.nan
        #smooth[-n:] = np.nan

    return smooth

def godin_shape():
    """
    Based on matlab code of 4/8/2013  Parker MacCready
    Returns a 71 element numpy array that is the weights
    for the Godin 24-24-25 tildal averaging filter. This is the shape given in
    Emery and Thomson (1997) Eqn. (5.10.37)
    ** use ONLY with hourly data! **
    """
    import numpy as np
    k = np.arange(12)
    filt = np.NaN * np.ones(71)
    filt[35:47] = (0.5/(24*24*25))*(1200-(12-k)*(13-k)-(12+k)*(13+k))
    k = np.arange(12,36)
    filt[47:71] = (0.5/(24*24*25))*(36-k)*(37-k)
    filt[:35] = filt[:35:-1]
    return filt

def hanning_shape(n=40):
    """
    Returns a Hanning window of the specified length.
    """
    import numpy as np
    ff = np.cos(np.linspace(-np.pi,np.pi,n+2))[1:-1]
    filt = (1 + ff)/2
    filt = filt / filt.sum()
    return filt

def ncd(fn_ds, pat=''):
    """
    Prints info on varibles in a NetCDF file or NetCDF dataset. Accepts a string 'pat' that
    can be used to filter the output.

    Example: zfun.ncd(fn_ds, pat='time')
    """
    # determine input type
    import netCDF4 as nc
    if isinstance(fn_ds, nc.Dataset):
        ds = fn_ds
    elif isinstance(fn_ds, str):
        try:
            ds = nc.Dataset(fn_ds)
        except:
            print('Input was not a NetCDF file')
            return
    else:
        print('Bad input type')
        return # exit the function

    # print information
    for vn in ds.variables:
        if len(pat) > 0:
            if pat in vn:
                print(ds.variables[vn])
        else:
            print(ds.variables[vn])

def auto_lims(fld):
    """
    A convenience function for automatically setting color limits.
    Input: a numpy array (masked OK)
    Output: tuple of good-guess colorsclae limits for a pcolormesh plot.
    """
    import numpy as np
    flo = np.floor(fld.mean() - fld.std())
    fhi = np.ceil(fld.mean() + fld.std())
    return (flo, fhi)

def earth_rad(lat_deg):
    """
    Calculate the Earth radius (m) at a latitude
    (from http://en.wikipedia.org/wiki/Earth_radius) for oblate spheroid

    INPUT: latitude in degrees

    OUTPUT: Earth radius (m) at that latitute
    """
    a = 6378.137 * 1000; # equatorial radius (m)
    b = 6356.7523 * 1000; # polar radius (m)
    import numpy as np
    cl = np.cos(np.pi*lat_deg/180)
    sl = np.sin(np.pi*lat_deg/180)
    RE = np.sqrt(((a*a*cl)**2 + (b*b*sl)**2) / ((a*cl)**2 + (b*sl)**2))
    return RE





