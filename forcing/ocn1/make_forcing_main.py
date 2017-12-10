# -*- coding: utf-8 -*-
"""
This is the main program for making the OCN forcing file.

It is designed to work with the new hycom1 archive of extracted files,
as well as regular forecasts.

2017.12.10 I added a planB flag to the forecast case
which copies the clm file from the previous day and
makes the final time a day later.

*******************************

To run from the command line in LiveOcean/driver/:
    
./driver_forcing1.sh -g cas1 -t base -f ocn1 -r backfill -0 20130101 -1 20130101

To test in python on mac:

cd ~/Documents/LiveOcean/forcing/ocn1

This is a normal case (no gaps)
run make_forcing_main.py -g cas1 -t base -r backfill -d 2013.01.01

This has a huge gap, but still works:
run make_forcing_main.py -g cas1 -t base -r backfill -d 2012.12.01

This will add the bio variables:
run make_forcing_main.py -g cas1 -t base -x lobio3 -r backfill -d 2013.01.01

Or try a forecast
run make_forcing_main.py -g cascadia1 -t base -r forecast

"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

#%% ****************** CASE-SPECIFIC CODE *****************

from datetime import datetime, timedelta
import shutil
import pickle
import netCDF4 as nc

import zrfun
import Ofun
import Ofun_nc
from importlib import reload
reload(Ofun)
reload(Ofun_nc)

start_time = datetime.now()

vnl_full = ['ssh','s3d','t3d','u3d','v3d']
exnum = '91.2'
testing = True

if testing == False:
    planB = False
elif testing == True:
    planB = True

if (Ldir['run_type'] == 'forecast') and (planB == False):    
    h_out_dir = out_dir = Ldir['LOogf_fd']      
    print('** START getting catalog')
    # create a list of url's of the preprocessed HYCOM files for this forecast
    try:
        fn_list = Ofun.get_hycom_file_list(exnum)
        print('** END getting catalog')
    except:
        planB = True
        
    if planB == False:
        # get a selection of the raw list (e.g. one per day)
        varf_dict, dt_list = Ofun.get_varf_dict(fn_list, Ldir)
        var_list = list(varf_dict.keys())    
        vnl_dict = {'ssh':['ssh'], 'ts3z':['s3d','t3d'], 'uv3z':['u3d','v3d']}
        #get the data and pack it in pickle files
        for vns in var_list:
            this_list = varf_dict[vns]
            if testing:
                this_list = [this_list[0]]
            for fn in this_list:
                a = Ofun.get_extraction(fn, vns)
                dts = datetime.strftime(a['dt'], '%Y.%m.%d')
                out_fn = h_out_dir + 'h' + dts + '.p'
                if os.path.exists(out_fn)== True:
                    print('  Opening ' + out_fn)
                    aa = pickle.load(open(out_fn, 'rb'))
                    for vn in vnl_dict[vns]:
                        aa[vn] = a[vn]
                        print('    adding ' + vn)
                else:
                    print('  Creating ' + out_fn)
                    aa = dict()
                    aa['dt'] = a['dt']
                    for vn in vnl_dict[vns]:
                        aa[vn] = a[vn]
                        print('    adding ' + vn)
                pickle.dump(aa, open(out_fn, 'wb'))
        h_in_dir = h_out_dir
        h_list0 = os.listdir(h_in_dir)
        h_list = [item for item in h_list0 if item[0] == 'h']
       
elif (Ldir['run_type'] == 'backfill') and (planB == False):
    # get a list of all available times
    h_in_dir = Ldir['data'] + 'hycom1/'        
    h_list0 = os.listdir(h_in_dir)
    h_list1 = [item for item in h_list0 if item[0] == 'h']
    # then find the index of the start of the current day
    # but if it is missing search for the most recent past one that exists
    keep_looking = True
    dt_now = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
    it0 = None
    counter = 0
    maxcount = 100 # this handles the biggest gap we have
    while keep_looking and counter < maxcount:
        dt_next = dt_now - timedelta(days=counter)
        dts_next = datetime.strftime(dt_next, '%Y.%m.%d')
        try:
            it0 = h_list1.index('h' + dts_next + '.p')
            keep_looking = False
            if counter > 0:
                print('Warning: Needed %d iterations' % (counter))
        except ValueError:
            counter += 1       
    # save the list of files
    if it0 == None:
        print('ERROR: no valid files found at nearby times')    
    else:        
        it_list = range(it0-2, it0+4)
        h_list = []
        for it in it_list:
            fn = h_list1[it]
            h_list.append(fn)
            #print(fn) # debugging
    # note that we don't actually write any of the files, because
    # they already exist in LiveOcean_data/hycom1/

if planB == False:
    # copy in the coordinates (assume those from hycom1 work)
    exnum1 = '91.2'
    c_in_dir = Ldir['data'] + 'hycom1/'
    c_out_dir = out_dir = Ldir['LOogf_fd']
    coords_dict = pickle.load(open(c_in_dir + 'coords_dict.p', 'rb'))
    coord_dict = dict()
    for vn in ['lon', 'lat', 'z']:
        coord_dict[vn] = coords_dict[exnum1][vn]
    pickle.dump(coord_dict, open(c_out_dir + 'coord_dict.p', 'wb'))    
        
    #%% filter in time
    fh_dir = Ldir['LOogf_fd']
    Ofun.time_filter(h_in_dir, h_list, fh_dir, Ldir)

    #%% extrapolate
    lon, lat, z, L, M, N, X, Y = Ofun.get_coords(fh_dir)
    a = os.listdir(fh_dir)
    aa = [item for item in a if item[:2]=='fh']
    for fn in aa:    
        print('-Extrapolating ' + fn)    
        in_fn = fh_dir + fn
        V = Ofun.get_extrapolated(in_fn, L, M, N, X, Y, z)
        pickle.dump(V, open(fh_dir + 'x' + fn, 'wb'))

    # and interpolate to ROMS format
    # get grid and S info
    G = zrfun.get_basic_info(Ldir['grid'] + 'grid.nc', only_G=True)
    S_info_dict = Lfun.csv_to_dict(Ldir['grid'] + 'S_COORDINATE_INFO.csv')
    S = zrfun.get_S(S_info_dict)
    # get list if files to work on
    a = os.listdir(fh_dir)
    aa = [item for item in a if item[:3]=='xfh']
    # HyCOM grid info
    lon, lat, z, L, M, N, X, Y = Ofun.get_coords(fh_dir)
    # load a dict of hycom fields
    dt_list = []
    count = 0
    c_dict = dict()
    for fn in aa:
        print('-Interpolating ' + fn + ' to ROMS grid')
        in_fn = fh_dir + fn
        b = pickle.load(open(in_fn, 'rb'))
        dt_list.append(b['dt'])
        c = Ofun.get_interpolated(G, S, b, lon, lat, z, N)   
        c_dict[count] = c
        count += 1
    
    #%% Write to ROMS forcing files
    nc_dir = Ldir['LOogf_f']
    Ofun_nc.make_clm_file(Ldir, nc_dir, fh_dir, c_dict, dt_list, S, G)
    
elif PlanB == True:
    ds_today = Ldir['date_string']
    dt_today = datetime.strptime(ds_today, format='%Y.%m.%d')
    dt_yesterday = dt_today - timedelta(days=1)
    ds_yesterday = datetime.strftime(dt_yesterday, format='%Y.%m.%d')
    clm_yesterday = (Ldir['LOog'] + 'f' + ds_yesterday + '/'
        + args.frc + '/' + 'ocean_clm.nc')
    clm_today = Ldir['LOogf_f'] + 'ocean_clm.nc'
    
    print(clm_yesterday)
    print(clm_today)
    # shutil.copyfile(clm_yesterday, clm_today)
    # ds = nc.Dataset(clm_today, 'a')

    
if False:
    # add bio variable if needed
    # if 'bio' in Ldir['ex_name']:
    #     Ofun_nc.add_bio(nc_dir)

    nc_dir = Ldir['LOogf_f']
    Ofun_nc.make_ini_file(nc_dir)
    Ofun_nc.make_bry_file(nc_dir)

    #%% prepare for finale
    import collections
    result_dict = collections.OrderedDict()
    time_format = '%Y.%m.%d %H:%M:%S'
    result_dict['start_time'] = start_time.strftime(time_format)
    end_time = datetime.now()
    result_dict['end_time'] = end_time.strftime(time_format)
    dt_sec = (end_time - start_time).seconds
    result_dict['total_seconds'] = str(dt_sec)

    ds = nc.Dataset(nc_dir + 'ocean_clm.nc')
    ot = ds['ocean_time'][:]
    ds.close()
    dt0 = Lfun.modtime_to_datetime(ot[0])
    dt1 = Lfun.modtime_to_datetime(ot[-1])

    result_dict['var_start_time'] = dt0.strftime(time_format)
    result_dict['var_end_time'] = dt1.strftime(time_format)

    nc_list = ['ocean_clm.nc', 'ocean_ini.nc', 'ocean_bry.nc']
    result_dict['result'] = 'success'
    for fn in nc_list:
        if os.path.isfile(nc_dir + fn):
            pass
        else:
           result_dict['result'] = 'fail'

    #%% ************** END CASE-SPECIFIC CODE *****************

    ffun.finale(result_dict, Ldir, Lfun)


