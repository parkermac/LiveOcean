"""
This is the main program for making the ATM forcing file.

Testing on my mac: need to use this day to find stored files:

run make_forcing_main.py -g cas6 -t v3 -d 2017.04.20

Note: we set rain to zero because its units are uncertain and
we don't currently use it in the simulations.
"""

import os; import sys
sys.path.append(os.path.abspath('../'))
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************
from datetime import datetime, timedelta
import time
import shutil
import netCDF4 as nc
import numpy as np
import seawater as sw
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
import matplotlib.path as mpath

import zfun
import zrfun

import atm_fun as afun
from importlib import reload
reload(afun)

start_time = datetime.now()

# where are files located, and other situational choices
do_d3 = True; do_d4 = True
if Ldir['lo_env'] == 'pm_mac':
    wrf_dir = Ldir['data'] + 'wrf/'
    Ldir['run_type'] == 'backfill'
    testing = False
    #do_d3 = True; do_d4 = False
    import matplotlib.pyplot as plt
    sys.path.append(os.path.abspath('../../plotting'))
    import pfun
else:
    wrf_dir = '/pmr2/darr/wrf_crons/wrfout/'
    testing = False
    
# creat list of hours
if Ldir['run_type'] == 'backfill':
    hr_vec = range(0,25)
elif Ldir['run_type'] == 'forecast':
    hr_max = Ldir['forecast_days'] * 24
    hr_vec = range(0, hr_max + 1)
    
# Create lists of input files.  These will be the full lists
# regardless of whether or not the files exist.
d_str = Ldir['date_string'].replace('.','')
in_dir = wrf_dir + d_str + '00/'
d2_list = []
d3_list = []
d4_list = []
forecast_hour = []
for hr in hr_vec:
    hr_str = ('0' + str(hr))[-2:]
    d2_list.append(in_dir + 'wrfout.ocean_d2.' + d_str + '00.f' + hr_str + '.0000')
    d3_list.append(in_dir + 'wrfout.ocean_d3.' + d_str + '00.f' + hr_str + '.0000')
    d4_list.append(in_dir + 'wrfout.ocean_d4.' + d_str + '00.f' + hr_str + '.0000')
    
# Create dicts that relate a filename to a time index
d2i_dict = {}
for i, v in enumerate(d2_list):
    d2i_dict[v] = i
        
# check for existence, and if any d2 are missing then exit
for fn in d2_list:
    if not os.path.isfile(fn):
        print('** Missing file: ' + fn)
        sys.exit() # this would be the place to invoke a Plan B
# for d3 and d4 just make sure we have the first one,
# so that we can get the grid
for fn in [d3_list[0]]:
    if not os.path.isfile(fn):
        print('** Missing file: ' + fn)
        do_d3 = False
for fn in [d4_list[0]]:
    if not os.path.isfile(fn):
        print('** Missing file: ' + fn)
        do_d4 = False

# create vector of time, in model format
dt0 = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
mod_time_list = []
for hr in hr_vec:
    dt = dt0 + timedelta(days=hr/24)
    mod_time = Lfun.datetime_to_modtime(dt)
    mod_time_list.append(mod_time)
mod_time_vec = np.array(mod_time_list)

# get model grid
gds = nc.Dataset(Ldir['grid'] + 'grid.nc')
lon = gds['lon_rho'][:]
lat = gds['lat_rho'][:]
gds.close()

# get WRF grid(s)
def get_wrf_grid(fn):
    wds = nc.Dataset(fn)
    lon = wds['XLONG'][:].squeeze()
    lat = wds['XLAT'][:].squeeze()
    if False:
        print('\n** ' + fn.split('/')[-1])
        vn_list = []
        for vn in wds.variables:
            vn_list.append(vn)
        print(vn_list)
    wds.close()
    # grid size info
    NR, NC = lon.shape
    jj = int(NR/2); ii = int(NC/2)
    dx_km, dd_deg = sw.dist(lat[jj,ii], [lon[jj,ii], lon[jj+1,ii+1]])
    return lon, lat, dx_km
# Note: lat, lon are only in the first file of the day (hour zero)
lon2, lat2, dx2_km = get_wrf_grid(d2_list[0])
if do_d3:
    try:
        lon3, lat3, dx3_km = get_wrf_grid(d3_list[0])
    except:
        do_d3 = False
if do_d4:
    try:
        # sometimes there are empty files
        lon4, lat4, dx4_km = get_wrf_grid(d4_list[0])
    except:
        do_d4 = False
    
# Limit varlist if testing
if testing == True:
    outvar_list = ['rain']
    #outvar_list = ['Pair','rain','swrad','lwrad_down','Tair','Qair']
else:
    outvar_list = afun.outvar_list
    
# initialize NetCDF output files, one for each variable
NR, NC = lon.shape
NT = len(mod_time_list)
ncformat = 'NETCDF3_64BIT_OFFSET'
tt0 = time.time()
nc_out_dict = {}
for vn in outvar_list:
    # name output file
    out_fn = Ldir['LOogf_f'] + vn + '.nc'
    nc_out_dict[vn] = out_fn
    # get rid of the old version, if it exists
    try:
        os.remove(out_fn)
    except OSError:
        pass # assume error was because the file did not exist
    foo = nc.Dataset(out_fn, 'w', format=ncformat)
    # create dimensions
    timename = afun.timename_dict[vn]
    foo.createDimension(timename, NT) # could use None
    foo.createDimension('eta_rho', NR)
    foo.createDimension('xi_rho', NC)
    # add time data
    vv = foo.createVariable(timename, float, (timename,))
    vv.units = 'seconds since 1970.01.01 UTC'
    vv[:] = mod_time_vec
    # add variable definition
    vv = foo.createVariable(vn, float, (timename, 'eta_rho', 'xi_rho'))
    vv.long_name = afun.longname_dict[vn]
    vv.units = afun.units_dict[vn]
    foo.close()
print('Initializing NetCDF files took %0.1f seconds' % (time.time() - tt0))
    

tt0 = time.time()
# find index to trim Eastern part of WRF fields
lon_max = lon[0,-1]
imax2 = zfun.find_nearest_ind(lon2[0,:], lon_max + .5)
lon2 = lon2[:,:imax2]; lat2 = lat2[:, :imax2]
if do_d3:
    imax3 = zfun.find_nearest_ind(lon3[0,:], lon_max + .5)
    lon3 = lon3[:,:imax3]; lat3 = lat3[:, :imax3]
if do_d4:
    imax4 = zfun.find_nearest_ind(lon4[0,:], lon_max + .5)
    lon4 = lon4[:,:imax4]; lat4 = lat4[:, :imax4]
# prepare coordinate arrays for interpolation
XY = np.array((lon.flatten(), lat.flatten())).T
XY2 = np.array((lon2.flatten(), lat2.flatten())).T
if do_d3:
    XY3 = np.array((lon3.flatten(), lat3.flatten())).T
if do_d4:
    XY4 = np.array((lon4.flatten(), lat4.flatten())).T
# find coordinate rotation matrices
def get_angle(lon, lat):
    NR, NC = lon.shape
    theta = np.nan * np.ones_like(lon)
    for jj in range(NR-1):
        junk, theta[jj,:-1] = sw.dist(lat[jj,:], lon[jj,:])
    theta[:,-1] = theta[:,-2]
    ca = np.cos(-np.pi*theta/180)
    sa = np.sin(-np.pi*theta/180)
    return ca, sa
ca2, sa2 = get_angle(lon2, lat2)
if do_d3:
    ca3, sa3 = get_angle(lon3, lat3)
if do_d4:
    ca4, sa4 = get_angle(lon4, lat4)
print('Manipulating grids took %0.1f seconds' % (time.time() - tt0))

# define regions for masking
def get_indices_in_polygon(plon_poly, plat_poly, lon, lat):
    # get Boolean mask array "M" that is true for points
    # in lon, lat that are in the polygon plon_poly, plat_poly
    V = np.ones((len(plon_poly),2))
    V[:,0] = plon_poly
    V[:,1] = plat_poly
    P = mpath.Path(V)
    M, L = lon.shape
    Rlon = lon.flatten()
    Rlat = lat.flatten()
    R = np.ones((len(Rlon),2))
    R[:,0] = Rlon
    R[:,1] = Rlat
    M = P.contains_points(R) # boolean
    M = M.reshape(lon.shape)
    return M
tt0 = time.time()
if do_d3:
    plon3_poly = np.concatenate((lon3[0,4:],lon3[:-5,-1],lon3[-5,4::-1],lon3[:-5:-1,4]))
    plat3_poly = np.concatenate((lat3[0,4:],lat3[:-5,-1],lat3[-5,4::-1],lat3[:-5:-1,4]))
    M3 = get_indices_in_polygon(plon3_poly, plat3_poly, lon, lat)
if do_d4:
    plon4_poly = np.concatenate((lon4[0,4:],lon4[:-5,-1],lon4[-5,4::-1],lon4[:-5:-1,4]))
    plat4_poly = np.concatenate((lat4[0,4:],lat4[:-5,-1],lat4[-5,4::-1],lat4[:-5:-1,4]))
    M4 = get_indices_in_polygon(plon4_poly, plat4_poly, lon, lat)
print('Make grid masks took %0.1f seconds' % (time.time() - tt0))

# get interpolation matrices
IM2 = cKDTree(XY2).query(XY); IM2 = IM2[1]
if do_d3:
    IM3 = cKDTree(XY3).query(XY); IM3 = IM3[1]
if do_d4:
    IM4 = cKDTree(XY4).query(XY); IM4 = IM4[1]

def gather_and_process_fields(fn, imax, ca, sa):
    # This is where we define any transformations to get from WRF to ROMS variables.
    ds = nc.Dataset(fn)
    iv_dict = dict()
    for ivn in afun.invar_list:
        # we trim fields to match the trimmed coordinate arrays
        iv_dict[ivn] = ds[ivn][0,:,:imax].squeeze()
    ds.close()
    # then convert to ROMS units/properties, still on the WRF grid
    # invar_list = ['Q2', 'T2', 'PSFC', 'U10', 'V10','RAINCV', 'RAINNCV', 'SWDOWN', 'GLW']
    # outvar_list = ['Pair','rain','swrad','lwrad_down','Tair','Qair','Uwind','Vwind']
    ov_dict = dict()
    for ovn in outvar_list:
        if ovn == 'Pair':
            # convert Pa to mbar
            ov_dict[ovn] = iv_dict['PSFC']/100 
        elif ovn == 'rain':
            # set this to zero because (a) we don't really understand the units
            # and (b) is it not used in the simulations at this point 2019.05.22
            ov_dict[ovn] = 0 * (iv_dict['RAINCV']+iv_dict['RAINNCV'])
        elif ovn == 'Tair':
            # convert K to C
            ov_dict[ovn] = iv_dict['T2'] - 273.15
        elif ovn == 'swrad':
            # account for reflection
            ov_dict[ovn] = iv_dict['SWDOWN'] * (1 - 0.1446)
        elif ovn == 'lwrad_down':
            # account for reflection
            ov_dict[ovn] = iv_dict['GLW']
        elif ovn == 'Qair':
            # calculate relative humidity [%]
             ov_dict[ovn] = afun.Z_wmo_RH(ov_dict['Pair'], ov_dict['Tair'], iv_dict['Q2'])
        elif ovn == 'Uwind':
            # % rotate velocity to E-W and N-S
            ov_dict[ovn] = ca*iv_dict['U10'] + sa*iv_dict['V10']
        elif ovn == 'Vwind':
            # % rotate velocity to E-W and N-S
            ov_dict[ovn] = ca*iv_dict['V10'] - sa*iv_dict['U10']
    return ov_dict
    
def interp_to_roms(ov_dict, outvar_list, XYn):
    # Interpolate to the ROMS grid, using nearest neighbor (about twice as fast as linear?)
    # Had we used linear interpolation it defaults to having nans outside the convex
    # hull of the data, but nearest neighbor interpolation fills everything, so instead we
    # use the masks M3 and M4 created above to decide where to add the data for finer grids.
    ovi_dict = dict()
    for ovn in outvar_list:
        v = ov_dict[ovn]
        ovi_dict[ovn] = griddata(XYn, v.flatten(), XY, method='nearest').reshape((NR,NC))
    return ovi_dict

def interp_to_roms_alt(ov_dict, outvar_list, IMn):
    ovi_dict = dict()
    for ovn in outvar_list:
        v = ov_dict[ovn].flatten()
        ovi_dict[ovn] = v[IMn].reshape((NR,NC))
    return ovi_dict

# MAIN TASK: loop over all hours

if testing == True:
    # 20 = about noon local time
    d2_list = d2_list[20:21]
    d3_list = d3_list[20:21]
    d4_list = d4_list[20:21]
dall_list = zip(d2_list, d3_list, d4_list)

for fn2, fn3, fn4 in dall_list:
    print('Working on ' + fn2.split('/')[-1] + ' and etc.')
    
    # if we are missing a d3 or d4 file then we don't do ANY of
    # that resolution after that
    if not os.path.isfile(fn3):
        print(' - missing ' + fn3)
        do_d3 = False
    if not os.path.isfile(fn4):
        print(' - missing ' + fn4)
        do_d4 = False
    
    tt0 = time.time()
    ov2_dict = gather_and_process_fields(fn2, imax2, ca2, sa2)
    ovi2_dict = interp_to_roms_alt(ov2_dict, outvar_list, IM2)
    print(' - d2: gather, process, and interp took %0.1f seconds' % (time.time() - tt0))
    
    if do_d3:
        try:
            tt0 = time.time()
            ov3_dict = gather_and_process_fields(fn3, imax3, ca3, sa3)
            ovi3_dict = interp_to_roms_alt(ov3_dict, outvar_list, IM3)
            print(' - d3: gather, process, and interp took %0.1f seconds' % (time.time() - tt0))
        except:
            print(' - could not process ' + fn3)
            do_d3 = False
    
    if do_d4:
        try:
            tt0 = time.time()
            ov4_dict = gather_and_process_fields(fn4, imax4, ca4, sa4)
            ovi4_dict = interp_to_roms_alt(ov4_dict, outvar_list, IM4)
            print(' - d4: gather, process, and interp took %0.1f seconds' % (time.time() - tt0))
        except:
            print(' - could not process ' + fn4)
            do_d4 = False
    
    tt0 = time.time()
    # combine the grids
    ovc_dict = dict()
    for ovn in outvar_list:
        v2 = ovi2_dict[ovn]
        v = v2.copy()
        if do_d3:
            v3 = ovi3_dict[ovn]
            v[M3] = v3[M3]
        if do_d4:
            v4 = ovi4_dict[ovn]
            v[M4] = v4[M4]
        if np.sum(np.isnan(v)) > 0:
            print('** WARNING Nans in combined output ' + ovn)
        ovc_dict[ovn] = v
    
    tt0 = time.time()
    # save to NetCDF
    tt = d2i_dict[fn2]
    for vn in outvar_list:
        fn = nc_out_dict[vn]
        foo = nc.Dataset(fn, 'a')
        foo[vn][tt,:,:] = ovc_dict[vn]
        foo.close()
    print(' - Write to NetCDF took %0.1f seconds' % (time.time() - tt0))
    
    
if testing == True:
    plt.close('all')
    lim_dict = dict(zip(afun.outvar_list, afun.lim_list))
    # plot some of the fields in the most recent ovi#_dicts
    for ovn in outvar_list:
        fig = plt.figure(figsize=(20,8))
        aa = [lon[0,0], lon[0,-1], lat[0,0], lat[-1,0]]
    
        ax = fig.add_subplot(131)
        ax.set_title('d2 ' + ovn)
        vmin, vmax = lim_dict[ovn]
        cs = plt.pcolormesh(lon, lat, ovi2_dict[ovn], cmap='rainbow', vmin=vmin, vmax=vmax)
        fig.colorbar(cs, ax=ax)
        pfun.dar(ax)
        pfun.add_coast(ax)
        ax.axis(aa)
    
        ax = fig.add_subplot(132)
        ax.set_title('Combined')
        fld = ovc_dict[ovn]
        cs = plt.pcolormesh(lon, lat, fld, cmap='rainbow', vmin=vmin, vmax=vmax)
        fig.colorbar(cs, ax=ax)
        pfun.dar(ax)
        pfun.add_coast(ax)
        ax.axis(aa)
    
        ax = fig.add_subplot(133)
        ax.set_title('Combined - d2')
        fld = ovc_dict[ovn] - ovi2_dict[ovn]
        vmax = np.max(np.abs([np.nanmax(fld),np.nanmin(fld)]))
        vmin = -vmax
        cs = plt.pcolormesh(lon, lat, fld, cmap='bwr', vmin=vmin, vmax=vmax)
        fig.colorbar(cs, ax=ax)
        pfun.dar(ax)
        pfun.add_coast(ax)
        ax.axis(aa)
    
        plt.show()
        

# ===========================================================================
# prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)

result_dict['result'] = 'success'
get_time = True
for vn in outvar_list:
    fn = nc_out_dict[vn]
    if os.path.isfile(fn):
        if get_time == True:
            ds = nc.Dataset(fn)
            mt0 = ds[afun.timename_dict[vn]][0]
            mt1 = ds[afun.timename_dict[vn]][-1]
            ds.close()
            dt0 = Lfun.modtime_to_datetime(float(mt0))
            dt1 = Lfun.modtime_to_datetime(float(mt1))
            result_dict['var_start_time'] = dt0.strftime(time_format)
            result_dict['var_end_time'] = dt1.strftime(time_format)
        get_time == False
    else:
       result_dict['result'] = 'fail'

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)