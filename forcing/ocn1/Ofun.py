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

hpth = os.path.abspath('../hycom1/')
if hpth not in sys.path:
    sys.path.append(hpth)
import hfun

import zfun
    
def get_hycom_file_list(exnum):
    """
    This parses the specified catalog.xml and returns a list containing
    the full paths to all the HYVOM NetCDF files in the rcurrent forecast.
    """   
    # look at Beautiful Soup?  Easier?
    import xml.etree.ElementTree as ET
    import urllib.request as U    
    # initiate the file list
    fn_list = []   
    # get the xml of the catalog
    xml_name = ('http://tds.hycom.org/thredds/catalog/GLBu0.08/expt_' + 
                exnum + '/forecasts/catalog.xml')
    try:
        xfile = U.urlopen(xml_name, timeout=30)
    except:
        print('problem getting xfile')
    tree = ET.parse(xfile)
    xfile.close()
    root = tree.getroot()
    rt = root.tag
    xmlns = rt[rt.find('{'): rt.find('}') + 1]
    # get the url prefix
    for e0 in root.findall('.//' + xmlns + 'service'):
        if e0.get('name') == 'ncdods':
            url_prefix = e0.get('base')   
    # get the remainder of the file paths and put them in a list
    for e0 in root.findall('.//' + xmlns + 'dataset'):
        if e0.get('urlPath') != None:
            fn_list.append(url_prefix + e0.get('urlPath'))            
    return fn_list
    
def get_varf_dict(fn_list, Ldir):  
    # Get list of daily filenames.
    # Names are like: .../hycom_glb_912_2017020400_t168_uv3z.nc 
    ssh_dict = dict()
    ts3z_dict = dict()
    uv3z_dict = dict()
    for fn in fn_list:
        fn1 = fn.split('/')[-1]
        fn2 = fn1.split('.')[0]
        parts = fn2.split('_')
        date_str = parts[-3]
        hour_str = parts[-2]
        var_str = parts[-1]
        hh = int(hour_str[1:])
        dt = datetime.strptime(date_str[:-2],'%Y%m%d') + timedelta(days=hh/24)
        if var_str == 'ssh':
            # by putting these in a dict we just retain the MOST RECENT
            # instance of a file for a given time, but they are NOT SORTED
            ssh_dict[dt] = fn
        elif var_str == 'ts3z':
            ts3z_dict[dt] = fn
        elif var_str == 'uv3z':
            uv3z_dict[dt] = fn
    # just save the times for which we have all three files
    # and which are at midnight
    dt_list2 = []        
    for dt in ssh_dict.keys():
        if (dt in ts3z_dict.keys()) and (dt in uv3z_dict.keys()) and (dt.hour==0):
            dt_list2.append(dt)
    #sort the datetime list            
    dt_list2.sort()
    # trim the list to get only what we need
    nd_f = Ldir['forecast_days']    
    dt0s = Ldir['date_string']
    dt0 = datetime.strptime(dt0s, '%Y.%m.%d')
    # assume we need two days before dt0, and two days after dt0+nd_f
    # but note that this ASSUMES we are using a 5 day window to filter in time.
    dt_list3 = []
    dt_low = dt0 - timedelta(days=2)
    dt_high = dt0 + timedelta(days=(nd_f+2))
    for dt in dt_list2:
        if dt>=dt_low and dt<=dt_high:
            dt_list3.append(dt)
    dt_list3.sort()    
    # pack results into varf_dict
    ssh_list = []
    ts3z_list = []
    uv3z_list = []
    for dt in dt_list3:
        ssh_list.append(ssh_dict[dt])
        ts3z_list.append(ts3z_dict[dt])
        uv3z_list.append(uv3z_dict[dt])
    varf_dict = dict()
    varf_dict['ssh'] = ssh_list    
    varf_dict['ts3z'] = ts3z_list    
    varf_dict['uv3z'] = uv3z_list    
    return (varf_dict, dt_list3)

def get_extraction(fn, var_name):
    # these are packed one time per file, so the time is encoded in "fn"    
    print('\n' + fn)
    # initialize an output dict
    out_dict = dict()    
    ds = nc.Dataset(fn)        
    # get the time in a meaningful format
    t = ds.variables['time'][:].squeeze() # expect shape (1,)
    # tu = ds.variables['time'].units
    # print(tu) # should be 'hours since 2000-01-01 00:00:00'
    t_origin = ds.variables['time'].time_origin
    from datetime import datetime
    from datetime import timedelta
    dt0 = datetime.strptime(t_origin, '%Y-%m-%d %H:%M:%S')       
    dt = (dt0 + timedelta(t/24.))   
    out_dict['dt'] = dt # datetime time of this snapshot    
    if var_name == 'ssh':
        out_dict['z'] = 0 # just for convenience
    else:
        # create z from the depth
        depth = ds.variables['depth'][:]
        z = -depth[::-1] # you reverse an axis with a -1 step!
        out_dict['z'] = z
        N = len(z)    
    # get full coordinates (vectors for the plaid grid)
    llon = ds.variables['lon'][:]
    llat = ds.variables['lat'][:]    
    # find indices of a sub region
    aa = hfun.get_extraction_limits()
    llon = llon - 360 # convert from 0:360 to -360:0 format
    i0 = zfun.find_nearest_ind(llon, aa[0])
    i1 = zfun.find_nearest_ind(llon, aa[1])
    j0 = zfun.find_nearest_ind(llat, aa[2])
    j1 = zfun.find_nearest_ind(llat, aa[3])   
    # and just get the desired region (these are vectors, not matrices)
    lon = llon[i0:i1]
    lat = llat[j0:j1]
    # and save them   
    out_dict['lon'] = lon
    out_dict['lat'] = lat    
    # start a timer
    tt0 = time.time()       
    # get the variables, renaming to be consistent with what we want
    # pack bottom to top
    # extrapolate horizontally to fill all space
    if var_name == 'ssh':
        ssh = ds.variables['surf_el'][0, j0:j1, i0:i1].squeeze()
        #print(str(ssh.shape))
        out_dict['ssh'] = ssh
    elif var_name == 'ts3z':
        t3d = ds.variables['water_temp'][0, 0:N, j0:j1, i0:i1].squeeze()
        t3d = t3d[::-1, :, :] # pack bottom to top
        out_dict['t3d'] = t3d
        s3d = ds.variables['salinity'][0, 0:N, j0:j1, i0:i1].squeeze()
        s3d = s3d[::-1, :, :]
        out_dict['s3d'] = s3d
    elif var_name == 'uv3z':
        u3d = ds.variables['water_u'][0, 0:N, j0:j1, i0:i1].squeeze()
        u3d = u3d[::-1, :, :]
        out_dict['u3d'] = u3d
        v3d = ds.variables['water_v'][0, 0:N, j0:j1, i0:i1].squeeze()
        v3d = v3d[::-1, :, :] # pack bottom to top
        out_dict['v3d'] = v3d
    print('  %0.2f sec to get %s' % ((time.time() - tt0), var_name))
    ds.close()
    return out_dict # now the keys of this dictionary are separate variables
                    
def time_filter(in_dir, h_list, out_dir, Ldir):
    """
    Filter the files that ended up in data_dir
    """
    print('\nFiltering in time\n')
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
        print('Using Hanning window')
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
            print(out_name)
            pickle.dump(aa, open(out_dir + out_name, 'wb'))
    else:
        print('Using block average')
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
        print(out_name0)
        pickle.dump(aa, open(out_dir + out_name0, 'wb'))
        # saving the last file
        dt1 = dt0 + timedelta(days=nd)
        dts1 = datetime.strftime(dt1, '%Y.%m.%d')
        out_name1 = 'fh' + dts1 + '.p'            
        aa['dt'] = dt1           
        print(out_name1)
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
    lonr = np.pi * lon / 180.0
    latr = np.pi * lat / 180.0
    # Create arrays of distance from the center (km) so that the
    # nearest neighbor extrapolation is based on physical distance
    RE = zfun.earth_rad(np.mean(lat))/1000 # radius of the Earth (km)
    mlatr = np.mean(latr)
    mlonr = np.mean(lonr)
    clat = np.cos(mlatr)
    X, Y = np.meshgrid(RE*clat*(lonr - mlonr), RE*(latr - mlatr))
    
    return (lon, lat, z, L, M, N, X, Y)

def checknan(fld):
    if np.isnan(fld).sum() > 0:
        print('WARNING: nans in data field')    

def extrap_nearest_to_masked(X, Y, fld, fld0=0):
    """
    INPUT: fld is a 2D array on spatial grid X, Y        
    OUTPUT: a numpy array of the same size with no mask
    and no missing values.        
    If input is a masked array:        
        * If is is ALL masked then return an array filled with fld0.         
        * If it is PARTLY masked use nearest neighbor interpolation to
        fill missing values, and then return data.        
        * If it is all unmasked then retun the data.    
    If input is not a masked array:        
        * Return the array.    
    """
    if np.ma.is_masked(fld):
        if fld.all() is np.ma.masked:
            print('  filling with ' + str(fld0))
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
    else:
        checknan(fld)
        return fld
            
def get_extrapolated(in_fn, L, M, N, X, Y, z):
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
    # some utility functions    

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
            v0 = v.min()
        elif vn == 's3d':
            v0 = v.max()   
        if vn in ['t3d', 's3d']:   
            for k in range(N):
                fld = v[k, :, :]
                fldf = extrap_nearest_to_masked(X, Y, fld, fld0=v0)
                V[vn][k, :, :] = fldf
        elif vn in ['u3d', 'v3d']:
            vv = v.copy()
            vv[v.mask] = 0
            V[vn] = vv.data            
    # Create ubar and vbar.
    # Note: this is slightly imperfect because the z levels are at the same
    # position as the velocity levels.
    dz = np.nan * np.ones((N, 1, 1))
    dz[1:, 0, 0]= np.diff(z)
    dz[0, 0, 0] = dz[1, 0, 0]
    dz3 = dz * np.ones_like(b['u3d']) # make dz a masked array
    b['ubar'] = np.sum(b['u3d']*dz3, axis=0) / np.sum(dz3, axis=0)
    b['vbar'] = np.sum(b['v3d']*dz3, axis=0) / np.sum(dz3, axis=0)
    for vn in ['ubar', 'vbar']:
        v = b[vn]
        vv = v.copy()
        vv[v.mask] = 0
        V[vn] = vv.data  
    # calculate potential temperature
    z3 = z.reshape((N,1,1))
    V['theta'] = seawater.ptmp(V['s3d'], V['t3d'], z3)
    
    return V

