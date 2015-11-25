"""
This is the main program for making the OCN forcing file.
"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************
from datetime import datetime

import Ofun
vnl_full = ['ssh','s3d','t3d','u3d','v3d']
# define the output location
nc_dir = Ldir['LOogf_fd']

if Ldir['run_type'] == 'forecast':
    print('** START getting catalog')
    # create a list of url's of the preprocessed HYCOM files for this forecast
    cc = 0
    cc_max = 10
    while cc < cc_max:
        try:
            print(' attempt number ' + str(cc))
            fn_list = Ofun.get_hycom_file_list()
            print('** DONE getting catalog')
            cc = cc_max
        except:
            cc += 1
    # get a selection of the raw list (e.g. one per day)
    var_list = ['ssh','ts3z','uv3z']    
    varf_dict = dict()
    for var_name in var_list:
        varf_dict[var_name] = Ofun.make_shortened_list(fn_list, var_name)
    # check that all the lists are the same length
    list_len = []
    for var_name in var_list:
        list_len.append(len(varf_dict[var_name]))
    if list_len.count(list_len[0]) != len(list_len):
        print('WARNING: Lists in varf_dict are different lengths!')
    else:
        NT = list_len[0] # the number of times, to use below   
    # use fewer times for testing
    if False:
        NT = 2       
    nt = 0
    try:
        while nt < NT:
            out_dict_allvar = dict()
            for var_name in var_list:
                fn = varf_dict[var_name][:][nt]           
                if var_name == 'ssh':
                    vnl = ['ssh']
                elif var_name == 'ts3z':
                    vnl = ['s3d','t3d']
                elif var_name == 'uv3z':
                    vnl = ['u3d','v3d']
                # this gets raw files (individual times)
                out_dict_multi = Ofun.get_extraction(fn, var_name)
                # but it is a bit messy because it often contains two variables           
                # so here we repack into a separate package for each variable
                out_dict = dict()
                # pack into a dict by variable name            
                for vn in vnl:                
                    items_list = ['z', 'dt', 'lon', 'lat', vn]
                    for item in items_list:
                        out_dict[item] = out_dict_multi[item]
                    # and add it to a dict of dicts    
                    out_dict_allvar[vn] = out_dict                      
            # create or append to the NetCDF files
            Ofun.fields_to_netcdf(out_dict_allvar, nt, nc_dir)
            nt += 1
        # do extrapolation and time filtering
        for vn in vnl_full:
            Ofun.process_extrap(vn, nc_dir)
            Ofun.process_time_filter(vn, nc_dir)
    except:
        # NOTE: 2015/06/15 the problem with this is that sometimes we get fields, but
        # not enough of them to do a three-day forecast.  In that case plan B
        # would have helped, but it is never executed.
        #
        # do this if geting the hycom forecast fields failed
        planB_flag = 0
        ndays = 1
        dt_now = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
        while (planB_flag == 0) and (ndays <= 6) :
            from datetime import timedelta
            dt_past = dt_now - timedelta(ndays)
            date_string_past = dt_past.strftime('%Y.%m.%d')
            Ldir['LOogfpast_f'] = (Ldir['LOo'] + Ldir['gtag'] +
                '/f' + date_string_past + '/' + Ldir['frc'] + '/')
            Ldir['LOogfpast_fd'] = Ldir['LOogfpast_f'] + 'Data/'            
            ssh_past = Ldir['LOogfpast_fd'] + 'ssh.nc'          
            import netCDF4 as nc    
            ds_past = nc.Dataset(ssh_past)
            t_past = ds_past.variables['dt'][:]
            ds_past.close()
            if t_past[-1] >= Lfun.datetime_to_modtime(dt_now):
                for vn in vnl_full:
                    fn_past = Ldir['LOogfpast_fd'] + vn + '.nc'
                    fn_new = Ldir['LOogf_fd'] + vn + '.nc'           
                    print('WARNING: unable to get todays forecast files!')
                    print('Instead copying:')
                    print('  ' + fn_past)
                    print('to:')
                    print('  ' + fn_new)
                    import shutil
                    shutil.copyfile(fn_past,fn_new)
                planB_flag = 1
            else:
                ndays += 1
                 
elif Ldir['run_type'] == 'backfill':
    import netCDF4 as nc
    import numpy as np
    from datetime import timedelta
    # set start and end times               
    dt0 = datetime.strptime(Ldir['date_string'],'%Y.%m.%d')
    dt1 = dt0 + timedelta(1)

    def make_smaller_netcdf(fn_in,fn_out, vn, nt0, nt1):
        import netCDF4 as nc
        ds = nc.Dataset(fn_in)
        lon = ds.variables['lon'][:]
        lat = ds.variables['lat'][:]
        tmod = ds.variables['tmod'][nt0:nt1+1]
        fld = ds.variables[vn + '_filt'][nt0:nt1+1]
        M, L = lon.shape   
        foo = nc.Dataset(nc_dir + vn + '.nc', 'w')
        foo.createDimension('t', None)               
        foo.createDimension('y', M)
        foo.createDimension('x', L)
        tmods = foo.createVariable('dt', float, ('t',))            
        lats = foo.createVariable('lat', float, ('y', 'x'))
        lons = foo.createVariable('lon', float, ('y', 'x'))
        lats[:] = lat
        lons[:] = lon
        tmods[:] = tmod
        if vn == 'ssh':
            flds = foo.createVariable(vn + '_filt', float, ('t', 'y', 'x'))
        else:
            z = ds.variables['z'][:]
            N = len(z)
            foo.createDimension('s', N)
            zs = foo.createVariable('z', float, ('s',))
            zs[:] = z
            flds = foo.createVariable(vn + '_filt', float, ('t', 's', 'y', 'x'))
        flds[:] = fld
        ds.close()
        foo.close()
    
    # find list of archived times (already extrapolated and filtered)
    vn_list = ['ssh', 't3d', 's3d', 'u3d', 'v3d']
    dir0 = Ldir['data'] + 'hycom_combined/'
    for vn in vn_list:
        fn = dir0 + vn + '.nc'
        ds = nc.Dataset(fn)
        tmod = ds.variables['tmod'][:]
        ds.close()
        dt_list = []
        for tt in tmod:
            dt_list.append(Lfun.modtime_to_datetime(tt))
        dt_arr = np.array(dt_list)
        dt00 = dt_arr[dt_arr <= dt0][-1]
        dt11 = dt_arr[dt_arr >= dt1][0]
        nt0 = dt_list.index(dt00)
        nt1 = dt_list.index(dt11)
        fn_out = nc_dir + vn + '.nc'
        print('Creating ' + fn_out)
        make_smaller_netcdf(fn,fn_out, vn, nt0, nt1)

# ************** END CASE-SPECIFIC CODE *****************

# run the code to create the forcing files
Lfun.run_worker(Ldir)

from datetime import datetime
print('MAIN end time = ' + str(datetime.now()))


