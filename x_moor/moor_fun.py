"""
Mooring extraction functions.
"""

import numpy as np
import netCDF4 as nc
import zfun

def get_itp_dict(sta_dict, G):
    # get interpolants
    itp_dict = dict()
    for sta_name in sta_dict.keys():
        # target position
        Lon = np.array(sta_dict[sta_name][0])
        Lat = np.array(sta_dict[sta_name][1])
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
        itp_dict[sta_name] = (Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
        
    return itp_dict

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
    
def start_netcdf(out_fn, N, NT, v0_list, v1_list, v2_list, v3_list_rho, v3_list_w, V_long_name, V_units):

    foo = nc.Dataset(out_fn, 'w')
    foo.createDimension('scalar', 1)
    foo.createDimension('s_rho', N)
    foo.createDimension('s_w', N+1)
    foo.createDimension('ocean_time', NT)

    for vv in v0_list:
        v_var = foo.createVariable(vv, float, ('scalar'))
        v_var.long_name = V_long_name[vv]
        v_var.units = V_units[vv]
    for vv in v1_list:
        v_var = foo.createVariable(vv, float, ('ocean_time',))
        v_var.long_name = V_long_name[vv]
        v_var.units = V_units[vv]
    for vv in v2_list:
        v_var = foo.createVariable(vv, float, ('ocean_time',))
        v_var.long_name = V_long_name[vv]
        v_var.units = V_units[vv]
    for vv in v3_list_rho:
        v_var = foo.createVariable(vv, float, ('ocean_time','s_rho'))
        v_var.long_name = V_long_name[vv]
        v_var.units = V_units[vv]
    for vv in v3_list_w:
        v_var = foo.createVariable(vv, float, ('ocean_time', 's_w'))
        v_var.long_name = V_long_name[vv]
        v_var.units = V_units[vv]

    foo.close()
    
def get_ij_good(lon, lat, mask, xvec, yvec, i0, j0):
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