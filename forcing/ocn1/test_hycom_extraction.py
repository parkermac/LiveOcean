#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 13:21:41 2017

@author: PM5

Stand-alone code to test the hycom extraction process for a daily forecast.
"""

from datetime import datetime, timedelta
import time

def get_hycom_file_list(exnum):
    """
    This parses the specified catalog.xml and returns a list containing
    the full paths to all the HYCOM NetCDF files in the current forecast.
    """   
    import xml.etree.ElementTree as ET
    from urllib.request import Request, urlopen
    from urllib.error import URLError
    from socket import timeout
    import time
             
    xml_name = ('http://tds.hycom.org/thredds/catalog/datasets/GLBu0.08/expt_' + 
                exnum + '/forecasts/catalog.xml')  
    req = Request(xml_name)
    counter = 1
    got_file = False
    while (counter <= 10) and (got_file == False):
        print('Attempting to get catalog XML, counter = ' + str(counter))    
        tt0 = time.time()
        try:
            xfile = urlopen(req, timeout=30)
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
        
    # initiate the file list
    fn_list = []           
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

def get_varf_dict(fn_list, date_string):  
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
    nd_f = 3    
    dt0s = date_string
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
    
    import netCDF4 as nc
    
    # these are packed one time per file, so the time is encoded in "fn"    
    print(fn)
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
    aa = get_extraction_limits()
    llon = llon - 360 # convert from 0:360 to -360:0 format
    i0 = find_nearest_ind(llon, aa[0])
    i1 = find_nearest_ind(llon, aa[1])
    j0 = find_nearest_ind(llat, aa[2])
    j1 = find_nearest_ind(llat, aa[3])   
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

def get_extraction_limits():
    # specify the sub region of hycom to extract
    aa = [-129, -121, 39, 51]
    return aa

def find_nearest_ind(array, value):
    import numpy as np
    # gives the index of the item in array that is closest to value
    idx = (np.abs(array-value)).argmin()
    return idx

# Here is the actual code we run, just focusing on seeing if we can
# successfully do the "get_extraction" step    
print('** START getting catalog')
date_string = datetime.now().strftime(format='%Y.%m.%d')
vnl_full = ['ssh','s3d','t3d','u3d','v3d']
exnum = '91.2'
testing = True
testing_shortlist = False
# create a list of url's of the preprocessed HYCOM files for this forecast
fn_list = get_hycom_file_list(exnum)
print('** END getting catalog')
# get a selection of the raw list (e.g. one per day)
varf_dict, dt_list = get_varf_dict(fn_list, date_string)
var_list = list(varf_dict.keys())    
vnl_dict = {'ssh':['ssh'], 'ts3z':['s3d','t3d'], 'uv3z':['u3d','v3d']}
#get the data and pack it in pickle files
for vns in var_list:
    this_list = varf_dict[vns]
    if testing_shortlist:
        this_list = [this_list[0]] # a shorter list
    for fn in this_list:
        if testing:
            print(fn)
        else:
            a = get_extraction(fn, vns)
            dts = datetime.strftime(a['dt'], '%Y.%m.%d')
            print('  Datetime of variable is ' + dts)


