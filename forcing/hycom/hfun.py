"""
Functions related getting archived HYCOM output.
"""

def get_hycom_past(fn, nt0, nt1):
 
    import netCDF4 as nc
           
    # initialize an output dict
    out_dict = dict()
    ds = nc.Dataset(fn)   
    # NOTE: time is handled by the calling function
    
    # create z from the depth
    depth = ds.variables['depth'][:]
    z = -depth[::-1] # you reverse an axis with a -1 step!
    N = len(z)
    out_dict['z'] = z
    
    # get full coordinates (vectors for the plaid grid)
    lon = ds.variables['lon'][:]
    lat = ds.variables['lat'][:]
    
    # find indices of a sub region
    import numpy as np
    lon = lon - 360 # convert to -360 to 0 format
    i0 = np.where(lon > -129)[0][0]
    i1 = np.where(lon < -121)[-1][-1]
    j0 = np.where(lat > 41)[0][0]
    j1 = np.where(lat < 51)[-1][-1]
    
    # and just get the desired region (these are vectors, not matrices)
    lon = lon[i0:i1]
    lat = lat[j0:j1]
    # and save them   
    out_dict['lon'] = lon
    out_dict['lat'] = lat
       
    import time
    tt0 = time.time() # start a timer
           
    # Get the variables, renaming to be consistent with what we want,
    # and re-pack bottom to top

    ssh = ds.variables['surf_el'][nt0:nt1+1, j0:j1, i0:i1].squeeze()
    out_dict['ssh'] = ssh

    t3d = ds.variables['water_temp'][nt0:nt1+1, 0:N, j0:j1, i0:i1].squeeze()
    t3d = t3d[:, ::-1, :, :] # pack bottom to top
    out_dict['t3d'] = t3d
    
    s3d = ds.variables['salinity'][nt0:nt1+1, 0:N, j0:j1, i0:i1].squeeze()
    s3d = s3d[:, ::-1, :, :]
    out_dict['s3d'] = s3d

    u3d = ds.variables['water_u'][nt0:nt1+1, 0:N, j0:j1, i0:i1].squeeze()
    u3d = u3d[:, ::-1, :, :]
    out_dict['u3d'] = u3d
    
    v3d = ds.variables['water_v'][nt0:nt1+1, 0:N, j0:j1, i0:i1].squeeze()
    v3d = v3d[:, ::-1, :, :] # pack bottom to top
    out_dict['v3d'] = v3d
              
    print '  ' + str(time.time() - tt0) + ' sec to get all variables'
            
    ds.close()    
    return out_dict
    
def hycom_dict_to_netcdf(out_dict, nc_dir, NT0, NT1):
    """
    Create or append the NetCDF files.
    
    NT0 and NT1 are the time index range in the NetCDF files where data from
    this out_dict will be placed.
    """
    import netCDF4 as nc
    import Lfun; reload(Lfun) # calling function must put this on the path
    
    vn_list = ['ssh', 't3d', 's3d', 'u3d', 'v3d']
    for vn in vn_list:        
        # loop over each variable, creating a separate NetCDF file for each    
        if NT0 == 0: # initialize if needed
            lon = out_dict['lon']
            lat = out_dict['lat']
            from numpy import meshgrid  
            Lon, Lat = meshgrid(lon,lat)
            M, L = Lon.shape   
            fld = out_dict[vn][:]
            # get rid of old versions
            import os
            try:
                os.remove(nc_dir + vn + '.nc')
            except OSError:
                pass # assume error was because the file did not exist 
            foo = nc.Dataset(nc_dir + vn + '.nc', 'w')
            foo.createDimension('t', None)               
            foo.createDimension('y', M)
            foo.createDimension('x', L)
            tmods = foo.createVariable('tmod', float, ('t',))            
            lats = foo.createVariable('lat', float, ('y', 'x'))
            lons = foo.createVariable('lon', float, ('y', 'x'))
            lats[:] = Lat
            lons[:] = Lon
            tmods[:] = out_dict['t'] # model time: seconds since 1/1/1970
            if vn == 'ssh':
                flds = foo.createVariable(vn, float, ('t', 'y', 'x'))
                flds[NT0:NT1+1,:,:] = fld
            else:
                z = out_dict['z']
                N = len(z)
                foo.createDimension('s', N)
                zs = foo.createVariable('z', float, ('s',))
                zs[:] = z
                flds = foo.createVariable(vn, float, ('t', 's', 'y', 'x'))
                flds[NT0:NT1+1,:,:,:] = fld
            foo.close()              
        else: # case when nt > 0, then we do not initialize
            foo = nc.Dataset(nc_dir + vn + '.nc', 'r+')
            fld = out_dict[vn][:]
            foo.variables['tmod'][NT0:NT1+1] = out_dict['t']
            if vn == 'ssh':
                foo.variables[vn][NT0:NT1+1,:,:] = fld
            else:
                foo.variables[vn][NT0:NT1+1,:,:,:] = fld                                 
            foo.close()  