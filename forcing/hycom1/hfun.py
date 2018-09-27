"""
Functions to use with the new HYCOM backfill code.

"""

# setup
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import zfun

import time
from datetime import datetime, timedelta

def get_extraction_limits():
    # specify the sub region of hycom to extract
    aa = [-129, -121, 39, 51]
    return aa

def get_dt_list(ds):
    # get the time in a meaningful format
    t_hycom = ds.variables['time'][:].squeeze()
    
    # hack
    import numpy as np
    t_hycom = np.concatenate((t_hycom, np.array([164208. + 48])))
    
    # tu = ds.variables['time'].units
    # print(' time units = ' + tu)
    # should be 'hours since 2000-01-01 00:00:00'
    t_origin = ds.variables['time'].time_origin
    dt00 = datetime.strptime(t_origin, '%Y-%m-%d %H:%M:%S')    
    # check the time reference
    if dt00 != datetime(2000,1,1,0,0,0):
        print('Warning: Unexpected time reference!')
        sys.stdout.flush()            
    dt_list = [] # make a list of datetime_values
    for tt in t_hycom:
        if tt/24 != int(tt/24):
            print('Warning: non-integer day!')
            sys.stdout.flush()
        dt_list.append(dt00 + timedelta(days=tt/24))
    return dt_list

def get_coordinates(ds):
    """
    Get lat, lon vectors, related indices, and z from the
    Dataset for a hycom archive.
    """
    # create z from the depth
    depth = ds.variables['depth'][:]
    # and pack bottom to top
    z = -depth[::-1]
    
    # get full coordinates (vectors for the plaid grid)
    llon = ds.variables['lon'][:]
    llat = ds.variables['lat'][:]
    
    # find indices of a sub region
    aa = get_extraction_limits()
    llon = llon - 360 # convert from 0:360 to -360:0 format
    i0 = zfun.find_nearest_ind(llon, aa[0])
    i1 = zfun.find_nearest_ind(llon, aa[1])
    j0 = zfun.find_nearest_ind(llat, aa[2])
    j1 = zfun.find_nearest_ind(llat, aa[3])
    
    # and just get the desired region (these are vectors, not matrices)
    lon = llon[i0:i1]
    lat = llat[j0:j1]
    
    coords = dict()
    coords['z'] = z
    coords['lon'] = lon
    coords['lat'] =lat
    coords['i0'] = i0
    coords['i1'] = i1
    coords['j0'] = j0
    coords['j1'] = j1
    
    return coords
    
def get_hycom_day(ds, nt, coords):
    """
    Extract the fields for a single day.
    """
    
    i0 = coords['i0']
    i1 = coords['i1']
    j0 = coords['j0']
    j1 = coords['j1']
    
    N = len(coords['z'])
           
    # initialize an output dict
    out_dict = dict()
               
    tt0 = time.time() # start a timer
           
    # Get the variables, renaming to be consistent with what we want,
    # and re-pack bottom to top

    ssh = ds.variables['surf_el'][nt, j0:j1, i0:i1].squeeze()    
    # test to see if we have good data at this time by looking
    # to see if ssh is all masked
    if ssh.mask.sum() == ssh.size:
        print('  No good data')
        sys.stdout.flush()
        out_dict['result'] = 'failure'
        return out_dict # and exit the function early
        
    else:
        # get and store all the fields
        out_dict['ssh'] = ssh
    
        t3d = ds.variables['water_temp'][nt, 0:N, j0:j1, i0:i1].squeeze()
        t3d = t3d[::-1, :, :]
        out_dict['t3d'] = t3d

        s3d = ds.variables['salinity'][nt, 0:N, j0:j1, i0:i1].squeeze()
        s3d = s3d[::-1, :, :]
        out_dict['s3d'] = s3d

        u3d = ds.variables['water_u'][nt, 0:N, j0:j1, i0:i1].squeeze()
        u3d = u3d[::-1, :, :]
        out_dict['u3d'] = u3d

        v3d = ds.variables['water_v'][nt, 0:N, j0:j1, i0:i1].squeeze()
        v3d = v3d[::-1, :, :]
        out_dict['v3d'] = v3d
        
        out_dict['result'] = 'success'              
        print('  %d sec to get all variables' % (int(time.time() - tt0)))
        sys.stdout.flush()          
        return out_dict