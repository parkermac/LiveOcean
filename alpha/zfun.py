"""
This module contains functions for extracting data from ROMS history files.

Parker MacCready
"""

import netCDF4 as nc
import numpy as np

def interp_scattered_on_plaid(x, y, xvec, yvec, u):
    """
    Gets values of the field u at locations (x,y).

    NOTE: this can also be used to interpolate to a plad grid.  Just pass it
    flattened versions of the new coordinate matrices, and then reshape
    the output.  Appears to be super fast

    All inputs and outputs are numpy arrays.

    Field u is defined on a plaid grid defined by vectors xvec and yvec.

    Locations (x,y) are defined by vectors whose elements specify
    individual points.  They are not a plaid grid, but can be
    scattered arbitrarily throughout the domain.

    Returns a vector ui the same length at x and y.

    Note that because it relies on "get_interpolant" out-of-bounds
    values default to nan.
    """
    # get interpolants
    xi0, xi1, xf = get_interpolant(x,xvec, extrap_nan=True)
    yi0, yi1, yf = get_interpolant(y,yvec, extrap_nan=True)

    # bi linear interpolation
    u00 = u[yi0,xi0]
    u10 = u[yi1,xi0]
    u01 = u[yi0,xi1]
    u11 = u[yi1,xi1]
    ui = (1-yf)*((1-xf)*u00 + xf*u01) + yf*((1-xf)*u10 + xf*u11)

    return ui

def get_interpolant(x, xvec, extrap_nan=False):
    """
    Returns info to allow fast interpolation.

    Input:
    x = data position(s) [1-D numpy array without nans]
    xvec = coordinate vector [1-D numpy array without nans]
        NOTE: xvec must be monotonically increasing

    *kwargs*
    Set extrap_nan=True to return nan for the fraction
    whenever an x value is outside of the range of xvec.
    The default is that
    if the x is out of the range of xvec it returns the
    interpolant for the first or last point.
    E.g. [0., 1., 0.] for x < xvec.min()

    Output: three 1-D numpy arrays of the same size as x
    i0 = index below [int]
    i1 = index above [int]
    fr = fraction [float]

    If the x is ON a point in xvec the default is to return
    the index of that point and the one above.

    """

    def itp_err(message='hi'):
        print('WARNING from get_interpolant(): ' + message)

    # input error checking
    if isinstance(x, np.ndarray) and isinstance(xvec, np.ndarray):
        pass # input type is good
    else:
        itp_err('Inputs must be numpy arrays')

    # some preconditioning of the input
    x = x.flatten()
    xvec = xvec.flatten()

    # more error checking
    if np.isnan(x).any():
        itp_err('nan found in x')
    if np.isnan(xvec).any():
        itp_err('nan found in xvec')
    if not np.all(np.diff(xvec) > 0):
        itp_err('xvec must be monotonic and increasing')

    nx = len(x)
    nxvec = len(xvec)

    X = x.reshape(nx, 1) # column vector
    xvec = xvec.reshape(1, nxvec)
    XVEC = xvec.repeat(nx, axis=0) # matrix

    # preallocate results arrays
    i0 = np.zeros(nx, dtype=int)
    i1 = np.zeros(nx, dtype=int)
    fr = np.zeros(nx, dtype=float)

    # calculate index columns
    mask = X >= XVEC
    # the above line broadcasts correctly even if nx = nxvec
    # because we forced X to be a column vector
    i0 = mask.sum(axis=1) - 1

    # these masks are used to handle values of x beyond the range of xvec
    lomask = i0 < 0
    himask = i0 > nxvec - 2
    i0[lomask] = 0
    i0[himask] = nxvec - 2
    i1 = i0 + 1

    # compute the fraction
    xvec0 = xvec[0,i0]
    xvec1 = xvec[0,i1]
    fr = (x - xvec0)/(xvec1 - xvec0)

    # fractions for out of range x
    if extrap_nan == False:
        fr[lomask] = 0.
        fr[himask] = 1.
    elif extrap_nan == True:
        fr[lomask] = np.nan
        fr[himask] = np.nan

    return i0, i1, fr

def find_nearest(array, value):
    # gives the itme in array that is closest to value
    idx = (np.abs(array-value)).argmin()
    return array[idx]

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

def get_S(S_info_dict):
    """
    Code to calculate S-coordinate vectors from the parameters
    in S_COORDINATE_INFO.csv.

    Need to check this carefully against the matlab version.

    PM 7/7/2016
    """

    #% recoded for python on 7/7/2016 from:
    #
    #% Z_scoord.m  5/21/2007  Parker MacCready
    #%
    #% this creates the structure S, which would be used for example by
    #% Z_s2z.m, given basic grid parameters
    #%
    #% IMPORTANT: note that this is the Song & Haidvogel 1994 stretching
    #% function, ROMS Vstretching = 1. If one chooses a different stretching
    #% function (as of March 2011 Vstretching can be 1, 2, or 3) then this code
    #% must be updated to include the proper stretching transformation!
    #%
    #% edited by DAS to include more things in S stucture
    #% edited by SNG March 2011 to include all of the current available ROMS
    #% stretching functions, 1-4 see:
    #% https://www.myroms.org/wiki/index.php/Vertical_S-coordinate#Vertical_Stretching_Functions

    S = dict()
    for item in S_info_dict.keys():
        if item in ['N', 'VSTRETCHING', 'VTRANSFORM']:
            S[item.title()] = int(S_info_dict[item])
        elif item in ['TCLINE', 'THETA_S', 'THETA_B']:
            S[item.lower()] = float(S_info_dict[item])
        else:
            pass

    N = S['N']
    Vstretching = S['Vstretching']
    Vtransform = S['Vtransform']
    tcline = S['tcline']
    theta_s = S['theta_s']
    theta_b = S['theta_b']

    hmin = 3 # a placeholder since I don't plan to use this.

    if Vtransform == 1:
        hc = min(hmin,tcline)
    elif Vtransform == 2:
        hc = tcline

    S['hc'] = hc

    sc_r = (np.linspace(-(N-1), 0, N) - 0.5)/N
    sc_w = np.linspace(-N, 0, N+1)/N

    S['sc_r'] = sc_r
    S['sc_w'] = sc_w

    if Vstretching == 1:
        if theta_s != 0:
            cff1 = 1/np.sinh(theta_s)
            cff2 = 0.5/np.tanh(0.5*theta_s)
            Cs_r = ( (1-theta_b)*cff1*np.sinh(theta_s*sc_r)
                    + theta_b*( cff2*np.tanh(theta_s*(sc_r + 0.5)) - 0.5 ) )
            Cs_w = ( (1-theta_b)*cff1*np.sinh(theta_s*sc_w)
                    + theta_b*( cff2*np.tanh(theta_s*(sc_w + 0.5)) - 0.5 ) )
        else:
            Cs_r = sc_r
            Cs_w = sc_w
    elif Vstretching == 2:
        alpha = 1
        beta = 1
        if theta_s!=0 and theta_b!=0:
            Csur = (1-np.cosh(theta_s*sc_r))/(np.cosh(theta_s)-1)
            Cbot = ((np.sinh(theta_b*(sc_r+1)))/(np.sinh(theta_b)))-1
            u = ((sc_r+1)**alpha)*(1+(alpha/beta)*(1-((sc_r+1)**beta)))
            Cs_r = u*Csur+(1-u)*Cbot
            Csur_w = (1-np.cosh(theta_s*sc_w))/(np.cosh(theta_s)-1)
            Cbot_w = ((np.sinh(theta_b*(sc_w+1)))/(np.sinh(theta_b)))-1
            u_w = ((sc_w+1)**alpha)*(1+(alpha/beta)*(1-((sc_w+1)**beta)))
            Cs_w = u_w*Csur_w+(1-u_w)*Cbot_w
        else:
            Cs_r = sc_r
            Cs_w = sc_w
    elif Vstretching == 3:
        # Geyer function for high bbl resolution in shallow applications
        gamma = 3
        Csur = -(np.log(np.cosh(gamma*abs(sc_r)**theta_s)))/np.log(np.cosh(gamma))
        Cbot = ((np.log(np.cosh(gamma*(sc_r+1)**theta_b)))/np.log(np.cosh(gamma)))-1
        mu = 0.5*(1-np.tanh(gamma*(sc_r+0.5)))
        Cs_r = mu*Cbot+(1-mu)*Csur
        Csur_w = -(np.log(np.cosh(gamma*abs(sc_w)**theta_s)))/np.log(np.cosh(gamma))
        Cbot_w = ((np.log(np.cosh(gamma*(sc_w+1)**theta_b)))/np.log(np.cosh(gamma)))-1
        mu_w = 0.5*(1-np.tanh(gamma*(sc_w+0.5)))
        Cs_w = mu_w*Cbot_w+(1-mu_w)*Csur_w
    elif Vstretching == 4:
        # newest ROMS default as of March 2011 (theta_s between 0 and 10,
        # theta_b between 0 and 4)
        if theta_s>0:
            Cs_r = (1-np.cosh(theta_s*sc_r))/(np.cosh(theta_s)-1)
            Cs_w = (1-np.cosh(theta_s*sc_w))/(np.cosh(theta_s)-1)
        elif theta_s<=0:
            Cs_r = -(sc_r**2)
            Cs_w = -(sc_w**2)
        if theta_b > 0:
            Cs_r = (np.exp(theta_b*Cs_r)-1)/(1-np.exp(-theta_b))
            Cs_w = (np.exp(theta_b*Cs_w)-1)/(1-np.exp(-theta_b))

    S['Cs_r'] = Cs_r
    S['Cs_w'] = Cs_w

    return S

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
    # input error checking (seems clumsy)
    if ( (type(h) != np.ndarray)
        or (type(zeta) not in [np.ndarray, np.ma.core.MaskedArray])
        or (type(S) != dict) ):
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
    ff = np.cos(np.linspace(-np.pi,np.pi,n+2))[1:-1]
    filt = (1 + ff)/2
    filt = filt / filt.sum()
    return filt

def ncd(fn_ds, pat=''):
    """
    Prints info on varibles in a NetCDF file or NetCDF dataset.
    Accepts a string 'pat' that
    can be used to filter the output.

    Example: zfun.ncd(fn_ds, pat='time')
    """
    # determine input type
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

def earth_rad(lat_deg):
    """
    Calculate the Earth radius (m) at a latitude
    (from http://en.wikipedia.org/wiki/Earth_radius) for oblate spheroid

    INPUT: latitude in degrees

    OUTPUT: Earth radius (m) at that latitute
    """
    a = 6378.137 * 1000; # equatorial radius (m)
    b = 6356.7523 * 1000; # polar radius (m)
    cl = np.cos(np.pi*lat_deg/180)
    sl = np.sin(np.pi*lat_deg/180)
    RE = np.sqrt(((a*a*cl)**2 + (b*b*sl)**2) / ((a*cl)**2 + (b*sl)**2))
    return RE





