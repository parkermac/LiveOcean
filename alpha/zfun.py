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
    # gives the item in array that is closest to value
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_nearest_ind(array, value):
    # gives the index of the item in array that is closest to value
    idx = (np.abs(array-value)).argmin()
    return idx

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
    Input: 1D numpy array data
    Output: Array of the same size, filtered with Hanning window of length n,
        padded with nan's
    If n=1 it just returns matrix you gave it
    """
    if n == 1:
        smooth = data
    else:
        filt = hanning_shape(n=n)
        n = np.ceil(len(filt)/2).astype(int)
        smooth = np.convolve(data, filt, mode = 'same')
        smooth[:n] = np.nan
        smooth[-n:] = np.nan
    return smooth

def filt_godin(data):
    """
    Input: 1D numpy array of HOURLY values
    Output: Array of the same size, filtered with 24-24-25 Godin filter,
        padded with nan's
    """
    filt = godin_shape()
    filt = filt / filt.sum()
    n = np.ceil(len(filt)/2).astype(int)
    smooth = np.convolve(data, filt, mode = 'same')
    smooth[:n] = np.nan
    smooth[-n:] = np.nan
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





