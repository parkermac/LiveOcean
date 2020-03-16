"""
Functions for the new ocean forcing.

"""
import os
import sys
import netCDF4 as nc
from datetime import datetime, timedelta
import numpy as np
import time
import pickle
from scipy.spatial import cKDTree
import seawater

import Ofun_CTD

hpth = os.path.abspath('../hycom2/')
if hpth not in sys.path:
    sys.path.append(hpth)
import hfun

import zfun
import zrfun

def get_data_oneday(this_dt, fn_out):
    """"
    Code to get hycom data using the new FMRC_best file. It gets only a single time,
    per the new guidance from Michael McDonald at HYCOM, 202.03.16.
    """

    import os
    import netCDF4 as nc

    from urllib.request import urlretrieve
    from urllib.error import URLError
    from socket import timeout
    import time
    from datetime import datetime, timedelta

    testing = False

    # get rid of the old version, if it exists
    try:
        os.remove(fn_out)
    except OSError:
        pass # assume error was because the file did not exist

    dstr = this_dt.strftime('%Y.%m.%d')
    print('Working on ' + dstr)
    
    # time string in HYCOM format
    dstr_hy = this_dt.strftime('%Y-%m-%d-T00:00:00Z')
    
    
    # specify spatial limits
    aa = hfun.aa
    north = aa[3]
    south = aa[2]
    west = aa[0] + 360
    east = aa[1] + 360

    if testing == True:
        var_list = 'surf_el'
    else:
        var_list = 'surf_el,water_temp,salinity,water_u,water_v'

    # create the request url (added new url 2019.12.06)
    url = ('https://ncss.hycom.org/thredds/ncss/GLBy0.08/expt_93.0/FMRC/GLBy0.08_930_FMRC_best.ncd'+
        '?var='+var_list +
        '&north='+str(north)+'&south='+str(south)+'&west='+str(west)+'&east='+str(east) +
        '&time_start='+dstr_hy+'&time_end='+dstr_hy +
        '&addLatLon=true&accept=netcdf4')
        
    print(url)

    # get the data and save as a netcdf file
    counter = 1
    got_file = False
    while (counter <= 3) and (got_file == False):
        print('Attempting to get data, counter = ' + str(counter))
        tt0 = time.time()
        try:
            (a,b) = urlretrieve(url,fn_out)
            # a is the output file name
            # b is a message you can see with b.as_string()
        except URLError as e:
            if hasattr(e, 'reason'):
                print(' *We failed to reach a server.')
                print(' -Reason: ', e.reason)
            elif hasattr(e, 'code'):
                print(' *The server couldn\'t fulfill the request.')
                print(' -Error code: ', e.code)
        except timeout:
            print(' *Socket timed out')
        else:
            got_file = True
            print(' Worked fine')
        print(' -took %0.1f seconds' % (time.time() - tt0))
        counter += 1

    # check results
    ds = nc.Dataset(fn_out)
    print('\nVariables:')
    for vn in ds.variables:
        print('- '+vn)
    # get time info for the forecast
    t = ds['time'][0]
    if isinstance(t, np.ma.MaskedArray):
        th = t.data
    else:
        th = t
    tu = ds['time'].units
    # e.g. 'hours since 2018-11-20 12:00:00.000 UTC'
    # Warning: Brittle code below!
    ymd = tu.split()[2]
    hmss = tu.split()[3]
    hms = hmss.split('.')[0]
    hycom_dt0 = datetime.strptime(ymd + ' ' + hms, '%Y-%m-%d %H:%M:%S')
    this_dt = hycom_dt0 + timedelta(days=(th/24))
    print('Target time = ' + dstr)
    print('Actual time = ' + this_dt.strftime('%Y-%m-%d-T00:00:00Z'))
    ds.close()

def get_data(this_dt, fn_out, nd_f):
    """"
    Code to get hycom data using the new FMRC_best file.
    """

    import os
    import netCDF4 as nc

    from urllib.request import urlretrieve
    from urllib.error import URLError
    from socket import timeout
    import time
    from datetime import datetime, timedelta

    testing = False

    # specify output file
    # get rid of the old version, if it exists
    try:
        os.remove(fn_out)
    except OSError:
        pass # assume error was because the file did not exist

    # specify time limits
    # get today's  date string in LO format
    
    #dstr00 = datetime.now().strftime('%Y.%m.%d')
    dstr00 = this_dt.strftime('%Y.%m.%d')
    
    print('Working on ' + dstr00)
    # get time limits for forecast
    dt00 = datetime.strptime(dstr00, '%Y.%m.%d')
    dt0 = dt00 - timedelta(days=2)
    if testing == True:
        dt1 = dt00 - timedelta(days=1)
    else:
        dt1 = dt00 + timedelta(days=int(nd_f) + 2)
    # put them in ncss format
    dstr0 = dt0.strftime('%Y-%m-%d-T00:00:00Z')
    dstr1 = dt1.strftime('%Y-%m-%d-T00:00:00Z')
    print('- dt0 = ' + dstr0)
    print('- dt1 = ' + dstr1)
    
    # specify spatial limits
    aa = hfun.aa
    #aa = hfun.get_extraction_limits()
    north = aa[3]
    south = aa[2]
    west = aa[0] + 360
    east = aa[1] + 360

    if testing == True:
        var_list = 'surf_el'
    else:
        var_list = 'surf_el,water_temp,salinity,water_u,water_v'

    # create the request url (added new url 2019.12.06)
    #url = ('http://ncss.hycom.org/thredds/ncss/grid/GLBy0.08/expt_93.0/data/forecasts/FMRC_best.ncd'+
    url = ('https://ncss.hycom.org/thredds/ncss/GLBy0.08/expt_93.0/FMRC/GLBy0.08_930_FMRC_best.ncd'+
        '?var='+var_list +
        '&north='+str(north)+'&south='+str(south)+'&west='+str(west)+'&east='+str(east) +
        '&time_start='+dstr0+'&time_end='+dstr1+'&timeStride=8' +
        '&addLatLon=true&accept=netcdf4')
        
    print(url)

    # get the data and save as a netcdf file
    counter = 1
    got_file = False
    while (counter <= 3) and (got_file == False):
        print('Attempting to get data, counter = ' + str(counter))
        tt0 = time.time()
        try:
            (a,b) = urlretrieve(url,fn_out)
            # a is the output file name
            # b is a message you can see with b.as_string()
        except URLError as e:
            if hasattr(e, 'reason'):
                print(' *We failed to reach a server.')
                print(' -Reason: ', e.reason)
            elif hasattr(e, 'code'):
                print(' *The server couldn\'t fulfill the request.')
                print(' -Error code: ', e.code)
        except timeout:
            print(' *Socket timed out')
        else:
            got_file = True
            print(' Worked fine')
        print(' -took %0.1f seconds' % (time.time() - tt0))
        counter += 1

    # check results
    ds = nc.Dataset(fn_out)
    print('\nVariables:')
    for vn in ds.variables:
        print('- '+vn)
    print('Times:')
    htime = ds['time'][:]
    # get time info for the forecast
    t = ds['time'][:]
    tu = ds['time'].units
    # e.g. 'hours since 2018-11-20 12:00:00.000 UTC'
    # Warning: Brittle code below!
    ymd = tu.split()[2]
    hmss = tu.split()[3]
    hms = hmss.split('.')[0]
    hycom_dt0 = datetime.strptime(ymd + ' ' + hms, '%Y-%m-%d %H:%M:%S')
    # make a list of datetimes that are in the forecast
    hycom_dt_list = [] # datetimes
    hycom_iit_list = [] # indices
    iit = 0
    for th in t: # assume t is sorted
        this_dt = hycom_dt0 + timedelta(th/24)
        hycom_dt_list.append(this_dt)
        print('- '+str(this_dt))
        hycom_iit_list.append(iit)
        iit += 1
    ds.close()

def get_hnc_short_list(this_dt, Ldir):
    # initial experiment list
    hy_list = list(hfun.hy_dict.keys())
    hy_list.sort()
    # only use part of hy_list, splitting based on the change of gridsize
    # between hy5 and hy6
    ihy_split = hy_list.index('hy6')
    if this_dt <= datetime(2018,12,6):
        hy_list = hy_list[:ihy_split]
    elif this_dt >= datetime(2018,12,7):
        hy_list = hy_list[ihy_split:]
    # create a list of daily HYCOM NetCDF files in the archive
    hy_in_dir = Ldir['data'] + 'hycom2/'
    hnc_list = []
    for hy in hy_list:
        in_dir = hy_in_dir + hy + '/'
        hnc_list0 = os.listdir(in_dir)
        hnc_list1 = [in_dir + item for item in hnc_list0 if '.nc' in item]
        hnc_list1.sort()
        hnc_list += hnc_list1
    # remove repeated days
    seen = []
    hnc_unique_list = []
    for hnc in hnc_list:
        dd = hnc.split('/')[-1]
        if dd not in seen:
            seen.append(dd)
            hnc_unique_list.append(hnc)
    # then find the index of the start of the current day
    # but if it is missing search for the most recent past one that exists
    keep_looking = True
    dt_now = this_dt#datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
    it0 = None
    counter = 0
    maxcount = 100 # this handles the biggest gap we have
    while keep_looking and counter < maxcount:
        dt_next = dt_now - timedelta(days=counter)
        dts_next = datetime.strftime(dt_next, '%Y.%m.%d')
        try:
            it0 = [i for i, s in enumerate(hnc_unique_list) if dts_next in s]
            it0 = it0[0]
            keep_looking = False
            if counter > 0:
                print('Warning: Needed %d iterations' % (counter))
        except (ValueError, IndexError):
            counter += 1
    # save the list of files
    if it0 == None:
        print('ERROR: no valid files found at nearby times')
    else:
        it_list = range(it0-2, it0+4)
        hnc_short_list = []
        for it in it_list:
            if it < 0:
                print('ERROR: Requested time is out of range of the list')
                sys.exit()
            fn = hnc_unique_list[it]
            hnc_short_list.append(fn)
            
    return hnc_short_list

def convert_extraction_oneday(fn):
    
    testing = False
    
    # initialize an output dict
    out_dict = dict()
    ds = nc.Dataset(fn)
    
    print('opening ' + fn)
    
    # get time info for the forecast
    t = ds['time'][0]
    if isinstance(t, np.ma.MaskedArray):
        th = t.data
    else:
        th = t
    tu = ds['time'].units
    # e.g. 'hours since 2018-11-20 12:00:00.000 UTC'
    # Warning: Brittle code below!
    ymd = tu.split()[2]
    hmss = tu.split()[3]
    hms = hmss.split('.')[0]
    hycom_dt0 = datetime.strptime(ymd + ' ' + hms, '%Y-%m-%d %H:%M:%S')
    this_dt = hycom_dt0 + timedelta(days=(th/24))
    out_dict['dt'] = this_dt # datetime time of this snapshot
    print('- for dt = ' + str(this_dt))
    
    if testing == False:
        # create z from the depth
        depth = ds.variables['depth'][:]
        z = -depth[::-1] # you reverse an axis with a -1 step!
        out_dict['z'] = z
        N = len(z)
        
    # get full coordinates (vectors for the plaid grid)
    lon = ds.variables['lon'][:] - 360  # convert from 0:360 to -360:0 format
    lat = ds.variables['lat'][:]    
    # and save them   
    out_dict['lon'] = lon
    out_dict['lat'] = lat
    
    if testing == True:
        var_list = ['surf_el']
    else:
        var_list = ['surf_el,water_temp,salinity,water_u,water_v']
        
    # get the dynamical variables
    tt0 = time.time()
    for var_name in var_list:
        # Get the variables, renaming to be consistent with what we want,
        # and pack bottom to top
        np.warnings.filterwarnings('ignore') # suppress warning related to nans in fields
        if var_name == 'surf_el':
            ssh = ds['surf_el'][0, :, :]
            out_dict['ssh'] = ssh
        elif var_name == 'water_temp':
            t3d = ds['water_temp'][0, :, :, :]
            t3d = t3d[::-1, :, :] # pack bottom to top
            out_dict['t3d'] = t3d
        elif var_name == 'salinity':
            s3d = ds['salinity'][0, :, :, :]
            s3d = s3d[::-1, :, :]
            out_dict['s3d'] = s3d
        elif var_name == 'water_u':
            u3d = ds['water_u'][0, :, :, :]
            u3d = u3d[::-1, :, :]
            out_dict['u3d'] = u3d
        elif var_name == 'water_v':
            v3d = ds['water_v'][0, :, :, :]
            v3d = v3d[::-1, :, :] # pack bottom to top
            out_dict['v3d'] = v3d
        # print('  %0.2f sec to get %s' % ((time.time() - tt0), var_name))
    ds.close()
    
    return out_dict # the keys of this dictionary are separate variables

def convert_extraction(fn, iit):
    # initialize an output dict
    out_dict = dict()
    ds = nc.Dataset(fn)
    
    print('opening ' + fn)
    print('- for iit = ' + str(iit))
    
    # get datetime
    t = ds['time'][iit]
    tu = ds['time'].units
    # e.g. 'hours since 2018-11-20 12:00:00.000 UTC'
    # Warning: Brittle code below!
    ymd = tu.split()[2]
    hmss = tu.split()[3]
    hms = hmss.split('.')[0]
    hycom_dt0 = datetime.strptime(ymd + ' ' + hms, '%Y-%m-%d %H:%M:%S')
    dt = hycom_dt0 + timedelta(days=t/24)
    out_dict['dt'] = dt # datetime time of this snapshot
    
    print('- for dt = ' + str(dt))
    
    # create z from the depth
    depth = ds.variables['depth'][:]
    z = -depth[::-1] # you reverse an axis with a -1 step!
    out_dict['z'] = z
    N = len(z)
        
    # get full coordinates (vectors for the plaid grid)
    lon = ds.variables['lon'][:] - 360  # convert from 0:360 to -360:0 format
    lat = ds.variables['lat'][:]    
    # and save them   
    out_dict['lon'] = lon
    out_dict['lat'] = lat
        
    # get the dynamical variables
    tt0 = time.time()
    for var_name in ['surf_el', 'water_temp', 'salinity', 'water_u', 'water_v']:
        # Get the variables, renaming to be consistent with what we want,
        # and pack bottom to top
        np.warnings.filterwarnings('ignore') # suppress warning related to nans in fields
        if var_name == 'surf_el':
            ssh = ds['surf_el'][iit, :, :]
            out_dict['ssh'] = ssh
        elif var_name == 'water_temp':
            t3d = ds['water_temp'][iit, :, :, :]
            t3d = t3d[::-1, :, :] # pack bottom to top
            out_dict['t3d'] = t3d
        elif var_name == 'salinity':
            s3d = ds['salinity'][iit, :, :, :]
            s3d = s3d[::-1, :, :]
            out_dict['s3d'] = s3d
        elif var_name == 'water_u':
            u3d = ds['water_u'][iit, :, :, :]
            u3d = u3d[::-1, :, :]
            out_dict['u3d'] = u3d
        elif var_name == 'water_v':
            v3d = ds['water_v'][iit, :, :, :]
            v3d = v3d[::-1, :, :] # pack bottom to top
            out_dict['v3d'] = v3d
        # print('  %0.2f sec to get %s' % ((time.time() - tt0), var_name))
    ds.close()
    
    return out_dict # the keys of this dictionary are separate variables

def time_filter(in_dir, h_list, out_dir, Ldir):
    """
    Filter the files that ended up in data_dir
    """
    print('-Filtering in time')
    vl = ['ssh', 'u3d', 'v3d', 't3d', 's3d']
    dts0 = Ldir['date_string']
    dt0 = datetime.strptime(dts0, '%Y.%m.%d')
    nh = len(h_list)    
    # test for gaps in h_list
    no_gaps = True
    dtsg = h_list[0].strip('h').strip('.p')
    dtg0 = datetime.strptime(dtsg, '%Y.%m.%d')
    for hh in h_list[1:]:
        dtsg = hh.strip('h').strip('.p')
        dtg1 = datetime.strptime(dtsg, '%Y.%m.%d')
        if (dtg1-dtg0).days != 1:
            no_gaps = False
            print('** HAS GAPS **')
            break
        else:
            dtg0 = dtg1    
    fac_list_H = [12, 4, 3, 4, 12]
    # inverse weighting factors for a Hanning window of length 5
    nfilt = len(fac_list_H)
    nd_f = Ldir['forecast_days']
    nhmin_f = nfilt + nd_f
    nhmin_b = nfilt + 1
    rtp = Ldir['run_type']   
    if ((nh==nhmin_b and rtp=='backfill') or (nh>=nhmin_f and rtp=='forecast')) and no_gaps:
        print('--Using Hanning window')
        fac_list = fac_list_H
        for nt in range(nh - 4):
            n_center = nt + 2
            aa = dict()
            for n in range(nfilt):
                nn = n + nt
                fn = in_dir + h_list[nn]
                a = pickle.load(open(fn, 'rb'))
                for v in vl:
                    if n == 0:
                        aa[v] = a[v]/fac_list[n]
                    else:
                        aa[v] = aa[v] + a[v]/fac_list[n]
            out_name = 'f' + h_list[n_center]
            dts = out_name.strip('fh').strip('.p')
            dt = datetime.strptime(dts, '%Y.%m.%d')
            aa['dt'] = dt
            print('   ' + out_name)
            pickle.dump(aa, open(out_dir + out_name, 'wb'))
    else:
        print('--Using block average')
        # make a simple average and use it for everything
        fac_list = list(nh * np.ones(nh))
        aa = dict()
        for n in range(nh):
            fn = in_dir + h_list[n]
            a = pickle.load(open(fn, 'rb'))
            for v in vl:
                if n == 0:
                    aa[v] = a[v]/fac_list[n]
                else:
                    aa[v] = aa[v] + a[v]/fac_list[n]        
        if rtp == 'backfill':
            nd = 1
        else:
            nd = 3
        # saving the first file
        out_name0 = 'fh' + dts0 + '.p' 
        aa['dt'] = dt0           
        print('   ' + out_name0)
        pickle.dump(aa, open(out_dir + out_name0, 'wb'))
        # saving the last file
        dt1 = dt0 + timedelta(days=nd)
        dts1 = datetime.strftime(dt1, '%Y.%m.%d')
        out_name1 = 'fh' + dts1 + '.p'            
        aa['dt'] = dt1           
        print('   ' + out_name1)
        pickle.dump(aa, open(out_dir + out_name1, 'wb'))

def get_coords(in_dir):
    # get coordinate fields and sizes
    coord_dict = pickle.load(open(in_dir + 'coord_dict.p', 'rb'))
    lon = coord_dict['lon']
    lat = coord_dict['lat']
    z = coord_dict['z']
    L = len(lon)
    M = len(lat)
    N = len(z)
    # Create arrays of distance from the center (m) so that the
    # nearest neighbor extrapolation is based on physical distance
    Lon, Lat = np.meshgrid(lon,lat)
    X, Y = zfun.ll2xy(Lon, Lat, lon.mean(), lat.mean())
    return (lon, lat, z, L, M, N, X, Y)

def checknan(fld):
    if np.isnan(fld).sum() > 0:
        print('WARNING: nans in data field')    

def extrap_nearest_to_masked(X, Y, fld, fld0=0):
    """
    INPUT: fld is a 2D array (np.ndarray or np.ma.MaskedArray) on spatial grid X, Y
    OUTPUT: a numpy array of the same size with no mask
    and no missing values.        
    If input is a masked array:        
        * If it is ALL masked then return an array filled with fld0.         
        * If it is PARTLY masked use nearest neighbor interpolation to
        fill missing values, and then return data.        
        * If it is all unmasked then return the data.    
    If input is not a masked array:        
        * Return the array.    
    """
    # first make sure nans are masked
    if np.ma.is_masked(fld) == False:
        fld = np.ma.masked_where(np.isnan(fld), fld)
        
    if fld.all() is np.ma.masked:
        #print('  filling with ' + str(fld0))
        fldf = fld0 * np.ones(fld.data.shape)
        fldd = fldf.data
        checknan(fldd)
        return fldd
    else:
        # do the extrapolation using nearest neighbor
        fldf = fld.copy() # initialize the "filled" field
        xyorig = np.array((X[~fld.mask],Y[~fld.mask])).T
        xynew = np.array((X[fld.mask],Y[fld.mask])).T
        a = cKDTree(xyorig).query(xynew)
        aa = a[1]
        fldf[fld.mask] = fld[~fld.mask][aa]
        fldd = fldf.data
        checknan(fldd)
        return fldd

def get_extrapolated(in_fn, L, M, N, X, Y, lon, lat, z, Ldir, add_CTD=False):
    b = pickle.load(open(in_fn, 'rb'))
    vn_list = list(b.keys())    
    # check that things are the expected shape
    def check_coords(shape_tuple, arr_shape):
        if arr_shape != shape_tuple:
            print('WARNING: array shape mismatch')
    for vn in vn_list:
        if vn == 'dt':
            pass
        elif vn == 'ssh':
            check_coords((M, L), b[vn].shape)
        else:
            check_coords((N, M, L), b[vn].shape)    
    # creat output array and add dt to it.
    vn_list.remove('dt')
    V = dict()
    for vn in vn_list:
        V[vn] = np.nan + np.ones(b[vn].shape)
    V['dt'] = b['dt']    
    # extrapolate ssh
    vn = 'ssh'
    v = b[vn]
    vv = extrap_nearest_to_masked(X, Y, v)
    V[vn] = vv
    vn_list.remove('ssh')    
    # extrapolate 3D fields
    for vn in vn_list:
        v = b[vn]
        if vn == 't3d':
            v0 = np.nanmin(v)
        elif vn == 's3d':
            v0 = np.nanmax(v)
        if vn in ['t3d', 's3d']:
            print(' -- extrapolating ' + vn)
            if add_CTD==False:
                for k in range(N):
                    fld = v[k, :, :]
                    fldf = extrap_nearest_to_masked(X, Y, fld, fld0=v0)
                    V[vn][k, :, :] = fldf
            elif add_CTD==True:
                print(vn + ' Adding CTD data before extrapolating')
                Cast_dict, sta_df = Ofun_CTD.get_casts(Ldir)
                for k in range(N):
                    fld = v[k, :, :]
                    zz = z[k]
                    xyorig, fldorig = Ofun_CTD.get_orig(Cast_dict, sta_df,
                        X, Y, fld, lon, lat, zz, vn)
                    fldf = Ofun_CTD.extrap_nearest_to_masked_CTD(X,Y,fld,
                        xyorig=xyorig,fldorig=fldorig,fld0=v0)
                    V[vn][k, :, :] = fldf
        elif vn in ['u3d', 'v3d']:
            print(' -- extrapolating ' + vn)
            vv = v.copy()
            vv = np.ma.masked_where(np.isnan(vv), vv)
            vv[vv.mask] = 0
            V[vn] = vv.data
    # Create ubar and vbar.
    # Note: this is slightly imperfect because the z levels are at the same
    # position as the velocity levels.
    dz = np.nan * np.ones((N, 1, 1))
    dz[1:, 0, 0]= np.diff(z)
    dz[0, 0, 0] = dz[1, 0, 0]
    
    # account for the fact that the new hycom fields do not show up masked
    u3d = np.ma.masked_where(np.isnan(b['u3d']),b['u3d'])
    v3d = np.ma.masked_where(np.isnan(b['v3d']),b['v3d'])
    dz3 = dz * np.ones_like(u3d) # make dz a masked array
    b['ubar'] = np.sum(u3d*dz3, axis=0) / np.sum(dz3, axis=0)
    b['vbar'] = np.sum(v3d*dz3, axis=0) / np.sum(dz3, axis=0)
    
    for vn in ['ubar', 'vbar']:
        v = b[vn]
        vv = v.copy()
        vv = np.ma.masked_where(np.isnan(vv), vv)
        vv[vv.mask] = 0
        V[vn] = vv.data  
    # calculate potential temperature
    press_db = -z.reshape((N,1,1))
    V['theta'] = seawater.ptmp(V['s3d'], V['t3d'], press_db)    
    return V

def get_xyr(G, vn):
    # ROMS lon, lat arrays
    if vn in ['ssh', 'theta', 's3d']:
        xr = G['lon_rho']
        yr = G['lat_rho']
    elif vn in ['ubar', 'u3d']:
        xr = G['lon_u']
        yr = G['lat_u']
    elif vn in ['vbar', 'v3d']:
        xr = G['lon_v']
        yr = G['lat_v']
    else:
        print('Unknown variable name for get_xyr: ' + vn)
    return xr, yr

def get_zr(G, S, vn):
    # get z on the ROMS grids
    h = G['h']
    if vn in ['theta', 's3d']:
        zr = zrfun.get_z(h, 0*h, S, only_rho=True)
    elif vn in ['u3d']:    
        xru, yru = get_xyr(G, 'ubar')
        hu = zfun.interp_scattered_on_plaid(G['lon_u'], G['lat_u'],
                    G['lon_rho'][0,:], G['lat_rho'][:,0], h, exnan=False)
        hu = np.reshape(hu, G['lon_u'].shape)
        zr = zrfun.get_z(hu, 0*hu, S, only_rho=True)    
    elif vn in ['v3d']:    
        hv = zfun.interp_scattered_on_plaid(G['lon_v'], G['lat_v'],
                    G['lon_rho'][0,:], G['lat_rho'][:,0], h, exnan=False)
        hv = np.reshape(hv, G['lon_v'].shape)
        zr = zrfun.get_z(hv, 0*hv, S, only_rho=True)
    else:
        print('Unknown variable name for get_zr: ' + vn)
    return zr

def get_interpolated(G, S, b, lon, lat, z, N):
    # OBSOLETE (and slow!)
    
    # start input dict
    c = dict()
    # 2D fields
    for vn in ['ssh', 'ubar', 'vbar']:
        tt0 = time.time()
        xr, yr = get_xyr(G, vn)
        Mr, Lr = xr.shape
        u = b[vn]
        uu = zfun.interp_scattered_on_plaid(xr, yr, lon, lat, u)
        uu = np.reshape(uu, xr.shape)
        if np.isnan(uu).any():
            print('Warning: nans in output array for ' + vn)
        else:
            c[vn] = uu
        print(' --' + vn + ' xy interpolation took %0.1f seconds' % (time.time() - tt0))
        
    # 3D fields
    for vn in ['theta', 's3d', 'u3d', 'v3d']:
        tt0 = time.time()
        xr, yr = get_xyr(G, vn)
        zr = get_zr(G, S, vn)
        Nr, Mr, Lr = zr.shape
        v = b[vn]
        # create an intermediate array, which is on the ROMS lon, lat grid
        # but has the HyCOM vertical grid
        # get interpolants
        xi0, xi1, xf = zfun.get_interpolant(xr, lon, extrap_nan=True)
        yi0, yi1, yf = zfun.get_interpolant(yr, lat, extrap_nan=True)
        # bi linear interpolation
        u00 = v[:,yi0,xi0]
        u10 = v[:,yi1,xi0]
        u01 = v[:,yi0,xi1]
        u11 = v[:,yi1,xi1]
        vi = (1-yf)*((1-xf)*u00 + xf*u01) + yf*((1-xf)*u10 + xf*u11)
        vi = vi.reshape((N, Mr, Lr))
        print(' --' + vn + ' xy interpolation took %0.1f seconds' % (time.time() - tt0))
        
        tt0 = time.time()
        # make interpolants to go from the HyCOM vertical grid to the
        # ROMS vertical grid
        I0, I1, FR = zfun.get_interpolant(zr, z, extrap_nan=True)
        zrs = zr.shape    
        vif = vi.flatten()
        LL = np.tile(np.arange(Lr), Mr*Nr)
        MM = np.tile(np.repeat(np.arange(Mr), Lr), Nr)
        f0 = LL + MM*Lr + I0*Mr*Lr
        f1 = LL + MM*Lr + I1*Mr*Lr
        vv = (1-FR)*vif[f0] + FR*vif[f1]
        vv = vv.reshape(zrs)
        if np.isnan(vv).any():
            print('Warning: nans in output array for ' + vn)
        else:
            c[vn] = vv
        print(' --' + vn + ' z interpolation took %0.1f seconds' % (time.time() - tt0))
    return c

def get_zinds(h, S, z):
    # Precalculate the array of indices to go from HYCOM z to ROMS z.
    # This just finds the vertical index in HYCOM z for each ROMS z_rho
    # value in the whole 3D array, with the index being the UPPER one of
    # the two HYCOM z indices that any ROMS z falls between.
    tt0 = time.time()
    zr = zrfun.get_z(h, 0*h, S, only_rho=True)
    zrf = zr.flatten()
    zinds = np.nan * np.ones_like(zrf)
    if isinstance(z, np.ma.MaskedArray):
        z = z.data
    for ii in range(len(z)-1):
        zlo = z[ii]; zhi = z[ii+1]
        mask = (zrf>zlo) & (zrf<=zhi)
        zinds[mask] = ii+1 # this is where the UPPER index is enforced
    zinds = zinds.astype(int)
    if isinstance(zinds, np.ma.MaskedArray):
        zinds = zinds.data
    print(' --create zinds array took %0.1f seconds' % (time.time() - tt0))
    return zinds

def get_interpolated_alt(G, S, b, lon, lat, z, N, zinds):
    # This does the horizontal and vertical interpolation to get from
    # extrapolated, filtered HYCOM fields to ROMS fields.
    #
    # We use fast nearest neighbor interpolation as much as possible.
    # Also we interpolate everything to the ROMS rho grid, and then crudely
    # interpolate to the u and v grids at the last moment.  Much simpler.
    
    # start input dict
    c = {}
    
    # precalculate useful arrays that are used for horizontal interpolation
    tt0 = time.time()
    if isinstance(lon, np.ma.MaskedArray):
        lon = lon.data
    if isinstance(lat, np.ma.MaskedArray):
        lat = lat.data
    Lon, Lat = np.meshgrid(lon,lat)
    XYin = np.array((Lon.flatten(), Lat.flatten())).T
    XYr = np.array((G['lon_rho'].flatten(), G['lat_rho'].flatten())).T
    h = G['h']
    IMr = cKDTree(XYin).query(XYr)[1]
    print(' --create IMr tree took %0.1f seconds' % (time.time() - tt0))
    
    # 2D fields
    tt0 = time.time()
    for vn in ['ssh', 'ubar', 'vbar']:
        vv = b[vn].flatten()[IMr].reshape(h.shape)
        if vn == 'ubar':
            vv = (vv[:,:-1] + vv[:,1:])/2
        elif vn == 'vbar':
            vv = (vv[:-1,:] + vv[1:,:])/2
        vvc = vv.copy()
        # always a good idea to make sure dict entries are not just pointers
        # to arrays that might be changed later, hent the .copy()
        c[vn] = vvc
        checknan(vvc)
    print(' --2d var xy interpolation took %0.1f seconds' % (time.time() - tt0))
        
    # 3D fields
    verbose = False
    # create intermediate arrays which are on the ROMS lon_rho, lat_rho grid
    # but have the HYCOM vertical grid (N layers)
    F = np.nan * np.ones(((N,) + h.shape))
    vi_dict = {}
    for vn in ['theta', 's3d', 'u3d', 'v3d']:
        tt0 = time.time()
        FF = F.copy()
        for nn in range(N):
            vin = b[vn][nn,:,:].flatten()
            FF[nn,:,:] = vin[IMr].reshape(h.shape)
        if verbose:
            print('==== FF for ' + vn)
            print(FF[:,10,10])
        checknan(FF)
        vi_dict[vn] = FF
    print(' --3d var xy interpolation took %0.1f seconds' % (time.time() - tt0))
    
    # do the vertical interpolation from HYCOM to ROMS z positions
    for vn in ['theta', 's3d', 'u3d', 'v3d']:
        tt0 = time.time()
        vi = vi_dict[vn]
        if verbose:
            print('==== vi for ' + vn)
            print(vi[:,10,10])
        hinds = np.indices((S['N'], G['M'], G['L']))
        vvf = vi[zinds, hinds[1].flatten(), hinds[2].flatten()]
        vv = vvf.reshape((S['N'], G['M'], G['L']))
        vvc = vv.copy()
        if verbose:
            print('==== vvc for ' + vn)
            print(vvc[:,10,10])
        if vn == 'u3d':
            vvc = (vvc[:,:,:-1] + vvc[:,:,1:])/2
        elif vn == 'v3d':
            vvc = (vvc[:,:-1,:] + vvc[:,1:,:])/2
        checknan(vvc)
        c[vn] = vvc
        print(' --' + vn + ' z interpolation took %0.1f seconds' % (time.time() - tt0))
    return c


