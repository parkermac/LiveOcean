"""
Module of functions for making ocean forcing.
"""

def get_hycom_file_list():
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
    xml_name = 'http://tds.hycom.org/thredds/catalog/GLBu0.08/expt_91.1/forecasts/catalog.xml'
    try:
        xfile = U.urlopen(xml_name, timeout=30)
    except:
        print('problem getting xfile')
    tree = ET.parse(xfile)
    xfile.close()
    root = tree.getroot()
    
    # NOTE: you find "xmlns" by looking at any instance of e0.tag (the part contained in {})
    # (just type e0.tag at the command line) or by looking at the "xmlns" attribute
    # in the first line of the XML listing when viewed in Chrome.
    # Or, here we automate the procedure (best)...
    rt = root.tag
    xmlns = rt[rt.find('{'): rt.find('}') + 1]
    
    # NOTE: in order to successfully find things in the XML you have to look at
    # it first, which you do by going to the url given by "xml_name" above
    # in Firefox or Chrome (Chrome is most complete)
    
    # .// selects all subelements, on all levels beneath the current element.
    # For example, .//egg selects all egg elements in the entire tree.
    
    # get the url prefix
    for e0 in root.findall('.//' + xmlns + 'service'):
        if e0.get('name') == 'ncdods':
            url_prefix = e0.get('base')
    
    # get the remainder of the file paths and put them in a list
    for e0 in root.findall('.//' + xmlns + 'dataset'):
        if e0.get('urlPath') != None:
            fn_list.append(url_prefix + e0.get('urlPath'))
            
    return fn_list
    
def make_shortened_list(fn_list, which_var):
    varf_list = []
    for fn in fn_list:
        if ('_' + which_var) in fn:
            varf_list.append(fn)
        
    # find unique days
    #
    # this year_list thing is to ensure that our search strategy works across
    # year boundaries.  It assumes that the date string is uniquely found
    # in a file name by searching for a string like '_2014' for example.
    from datetime import datetime   
    this_year = int(datetime.now().strftime('%Y'))
    year_list = range(2000,this_year+2) # a list of integer years to try
    dd_prev = 'junk'
    dlist = []
    for fn in varf_list:
        for year in year_list:
            ds_try = fn.find('_' + str(year))
            if ds_try != -1:
                ds0 = ds_try + 1
                break
        ds1 = ds0 + 8
        dd = fn[ds0:ds1]
        if dd == dd_prev:
            pass
        else:
            dlist.append(dd)
        dd_prev = dd
        
    dlist0 = dlist[:-1] # all but the last
    dlast = dlist[-1] # the last, which goes into the future
    
    varf_list_short = []
    # get the 000 hour forecast for all past days    
    for ddd in dlist0:
        for fn in varf_list:
            if ddd in fn and '_t000' in fn:
                varf_list_short.append(fn)      
    # get daily values up to 7 days into the future
    hour_list = range(0,168+24,24) 
    for fn in varf_list:
        if dlast in fn:
            for hr in hour_list:
                if '_t' + ('000' + str(hr))[-3:] in fn:
                    varf_list_short.append(fn)
                    
    # return a list of url's for this variable, one file per day            
    return varf_list_short

def get_extraction(fn, var_name):
    # these are packed one time per file, so the time is encoded in "fn"
    import time 
    import numpy as np
    
    print('working on ' + fn)
    # initialize an output dict
    out_dict = dict()
    
    import netCDF4 as nc    
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
    lon = ds.variables['lon'][:]
    lat = ds.variables['lat'][:]
    
    # find indices of a sub region
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
    
    # start a timer
    tt0 = time.time()       
    # get the variables, renaming to be consistent with what we want
    # pack bottom to top
    # extrapolate horizontally to fill all space
    if var_name == 'ssh':
        ssh = ds.variables['surf_el'][0, j0:j1, i0:i1].squeeze()
        print(str(ssh.shape))
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
    print(str(time.time() - tt0) + ' sec to get ' + var_name)
    ds.close()
    return out_dict # now the keys of this dictionary are separate variables
    
def fields_to_netcdf(out_dict_allvar, nt, nc_dir):
    """
    Create or append the NetCDF files.
    """
    import netCDF4 as nc
    import Lfun
    # we are assuming Lfun is on our path because the path was added in the
    # calling function - is this a good idea?
    for vn in out_dict_allvar.keys():
        
        # also, we pass in one snapshot at a time, and this code concatenates
        # them (initializes if nt = 0)
        
        # get the results for this collection of variables
        out_dict = out_dict_allvar[vn]
        
        # time, in model units
        dt = out_dict['dt'] 
        t = Lfun.datetime_to_modtime(dt) 
        
        if vn == 'ssh':
            do_z = False
        else:
            do_z = True
        
        # loop over each variable, creating a separate NetCDF file for each    

        if nt == 0: # initialize if needed
            lon = out_dict['lon']
            lat = out_dict['lat']
            from numpy import meshgrid  
            Lon, Lat = meshgrid(lon,lat)
            M, L = Lon.shape   
            fld = out_dict[vn][:]
            # get rid of old versions (could use os.join for portability)
            # test this!!! may not be needed
            import os
            try:
                os.remove(nc_dir + vn + '.nc')
            except OSError:
                pass # assume error was because the file did not exist 
            foo = nc.Dataset(nc_dir + vn + '.nc', 'w')
            foo.createDimension('t', None)               
            foo.createDimension('y', M)
            foo.createDimension('x', L)
            dts = foo.createVariable('dt', float, ('t',))            
            lats = foo.createVariable('lat', float, ('y', 'x'))
            lons = foo.createVariable('lon', float, ('y', 'x'))
            lats[:] = Lat
            lons[:] = Lon
            dts[nt] = t
            if do_z:
                z = out_dict['z']
                N = len(z)
                foo.createDimension('s', N)
                zs = foo.createVariable('z', float, ('s',))
                zs[:] = z
                flds = foo.createVariable(vn, float, ('t', 's', 'y', 'x'))
                flds[nt,:,:,:] = fld
            else:
                flds = foo.createVariable(vn, float, ('t', 'y', 'x'))
                flds[nt,:,:] = fld
        else: # case when nt > 0, then we do not initialize
            fld = out_dict[vn][:]
            foo = nc.Dataset(nc_dir + vn + '.nc', 'r+')
            foo.variables['dt'][nt] = t
            if do_z:
                foo.variables[vn][nt,:,:,:] = fld
            else:
                foo.variables[vn][nt,:,:] = fld                   
        foo.close()

def process_extrap(vn, nc_dir):
    """
    Add an extrapolated version of the variable to the NetCDf file.
    """
    import sys
    import netCDF4 as nc
    import time
    foo = nc.Dataset(nc_dir + vn + '.nc', 'r+')
    fld = foo.variables[vn][:]
    fld_new = fld.copy()
    shp = fld.shape
    NT = shp[0]
    #NT = 1 # debugging
    if vn == 'ssh':
        NZ = 1
    else:
        NZ = shp[1]
    lon = foo.variables['lon'][0,:]
    lat = foo.variables['lat'][:,0]    
    # setup: create x, y
    tt1 = time.time()
    import numpy as np
    lonr = np.pi * lon / 180.0
    latr = np.pi * lat / 180.0
    # and create arrays of distance from the center (km) so that the
    # nearest neighbor extrapolation is based on physical distance
    RE = 6371.0 # radius of the Earth (km)
    mlatr = np.mean(latr)
    mlonr = np.mean(lonr)
    clat = np.cos(mlatr)
    x, y = np.meshgrid(RE*clat*(lonr - mlonr), RE*(latr - mlatr))
    print('  -- x,y setup %.3f seconds' % (time.time() - tt1))
    if vn == 'ssh':
        try:
            fve = foo.createVariable(vn + '_extrap', float, ('t', 'y', 'x'))
            fve[:] = fld_new
        except:
            pass # assume variable exists
    else:
        try:
            fve = foo.createVariable(vn + '_extrap', float, ('t', 's', 'y', 'x'))
            fve[:] = fld_new
        except:
            pass # assume variable exists

    for tt in range(NT):
        tt0 = time.time()
        print('  ' + vn + ':  extrapolating ' + str(tt) + ' of ' + str(NT))
        # new version: much faster, AFTER I moved the x,y setup OUT of the
        # time loop
        this_fld = fld_new[tt].squeeze()
        
        # trap to fix bad salinity fields
        if vn == 's3d':
            test_cell = this_fld[-1,40,40]
        if np.ma.count_masked(test_cell) == 1:
            print('fixing a bad salinity field')
            this_fld_prev = fld_new[tt-1].squeeze()
            this_fld_next = fld_new[tt+1].squeeze()
            this_fld = (this_fld_prev + this_fld_next)/2.0

        this_flde = horizontal_extrapolation(x, y, NZ, this_fld, vn) 
        print('   - took %.3f seconds' % (time.time() - tt0))
        sys.stdout.flush()
        foo.variables[vn + '_extrap'][tt] = this_flde
    print(' take a deep breath...')
    sys.stdout.flush()
    foo.close()
            
def horizontal_extrapolation(x, y, N, fld, fld_name):
    import numpy as np
    # do the extrapolation    
    from scipy.spatial import cKDTree 
    if fld_name == 'ssh':
        fldf = fld.copy() # initialize the "filled" field
        xygood = np.array((x[~fld.mask],y[~fld.mask])).T
        xybad = np.array((x[fld.mask],y[fld.mask])).T
        fldf[fld.mask] = fld[~fld.mask][cKDTree(xygood).query(xybad)[1]]
    elif fld_name in ['u3d', 'v3d']:
        fldf = fld.copy() 
        fldf[fldf.mask] = 0.0 # fill all missing values with 0.0 velocity           
    elif fld_name in ['t3d', 's3d']:
        fldf = fld.copy() # initialize  
        # extrapolate horizontally
        iz = 0
        while iz < N:            
            # pull out a slice and extrapolate
            fldi = fld[iz, :, :].squeeze()
            xygood = np.array((x[~fldi.mask],y[~fldi.mask])).T
            xybad = np.array((x[fldi.mask],y[fldi.mask])).T
            fldfi = fldi.copy()
            if len(xygood) > 0:
                fldfi[fldi.mask] = fldi[~fldi.mask][cKDTree(xygood).query(xybad)[1]]
            else: # what we do when there is no data in a layer
                if fld_name == 't3d':                   
                    fldfi = np.nanmin(fld)
                elif fld_name == 's3d':
                    fldfi = np.nanmax(fld)
            # save the extrapolated slice            
            fldf[iz, :, :] = fldfi
            iz += 1
    # return the results
    return fldf
                    
def process_time_filter(vn, nc_dir):
    """
    Add a filtered version of the variable to the NetCDf file.
    """
    import netCDF4 as nc
    foo = nc.Dataset(nc_dir + vn + '.nc', 'r+')
    fld = foo.variables[vn + '_extrap'][:]
    flds = fld.copy()
    sz = fld.shape
    NT = sz[0]
    if vn == 'ssh':
        try:
            fve = foo.createVariable(vn + '_filt', float, ('t', 'y', 'x'))
            fve[:] = flds
        except:
            pass # assume variable exists
    else:
        try:
            fve = foo.createVariable(vn + '_filt', float, ('t', 's', 'y', 'x'))
            fve[:] = flds
        except:
            pass # assume variable exists
    # make a smoothed version of fld, to remove inertial oscillations

    if NT > 5:
        if vn == 'ssh':    
            flds[2:-2,:,:] = (fld[:-4,:,:]/12. + fld[1:-3,:,:]/4. + fld[2:-2,:,:]/3.
            + fld[3:-1,:,:]/4. + fld[4:,:,:]/12.)
            flds[:2,:,:] = flds[2,:,:]
            flds[-2:,:,:] = flds[-3,:,:]
        else:
            flds[2:-2,:,:,:] = (fld[:-4,:,:,:]/12. + fld[1:-3,:,:,:]/4.
            + fld[2:-2,:,:,:]/3. + fld[3:-1,:,:,:]/4. + fld[4:,:,:,:]/12.)
            flds[:2,:,:,:] = flds[2,:,:,:]
            flds[-2:,:,:,:] = flds[-3,:,:,:]
    else:
        import numpy as np
        ff = np.mean(flds, axis=0, keepdims=True)
        flds[:] = ff[:] # broadcasting, I hope!
    
    foo.variables[vn + '_filt'][:] = flds
    foo.close()
