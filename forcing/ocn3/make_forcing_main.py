# -*- coding: utf-8 -*-
"""
This is the main program for making the OCN forcing file.

It is designed to work with the new hycom1 archive of extracted files,
as well as regular forecasts.

2017.12.10 I added a planB flag to the forecast case
which copies the clm file from the previous day and
makes the final time a day later.

2018.05.19 I added the add_CTD flag (and the Ofun_CTD module) to add CTD data
on a specified day only.

*******************************

To run from the command line in LiveOcean/driver/:
    
./driver_forcing2.sh -g cas4 -t v2 -f ocn3 -r backfill -0 20170101 -1 20170101

To test in python on mac:

# standard backfill
run make_forcing_main.py -g cas4 -t v2 -r backfill -d 2017.01.10

# backfill with Salish and coastal estuary IC's from CTD and other info
run make_forcing_main.py -g cas4 -t v2 -r backfill -d 2017.01.01

# today's forecast
run make_forcing_main.py -g cas4 -t v2 -r forecast

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
import numpy as np

import zrfun
import Ofun
import Ofun_nc
import Ofun_CTD
import Ofun_bio
from importlib import reload
reload(Ofun)
reload(Ofun_nc)
reload(Ofun_CTD)
reload(Ofun_bio)

start_time = datetime.now()

# defaults
testing = False
planB = False
add_CTD = False
do_bio = True
fix_NSoG = False

# *** automate when to set add_CTD to True ***
this_dt = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
if this_dt == datetime(2017,1,1):
    print('WARNING: adding CTD data to extrapolation!!')
    add_CTD = True
    
# *** automate fixing of NSoG ***
if Ldir['run_type'] == 'backfill':
    fix_NSoG = True

if (Ldir['run_type'] == 'forecast'):
    # this either gets new hycom files, or sets planB to True
    
    h_out_dir = Ldir['LOogf_fd']
    Lfun.make_dir(h_out_dir, clean=True)
    
    try:
        data_fn_out =  h_out_dir + 'data.nc'
        nd_f = np.ceil(Ldir['forecast_days'])
        Ofun.get_data(this_dt, data_fn_out, nd_f)
        
        ds = nc.Dataset(data_fn_out)
        NT = len(ds['time'][:])
        ds.close()
        
        for iit in range(NT):
            a = dict()
            # get a dict of extractions at this time
            a = Ofun.get_extraction_new(data_fn_out, iit)
            dts = datetime.strftime(a['dt'], '%Y.%m.%d')
            out_fn = h_out_dir + 'h' + dts + '.p'
            pickle.dump(a, open(out_fn, 'wb'))

        h_in_dir = h_out_dir
        h_list0 = os.listdir(h_in_dir)
        h_list0.sort()
        h_list = [item for item in h_list0 if item[0] == 'h']
        
    except:
        print('- error getting forecast files')
        planB = True
       
elif (Ldir['run_type'] == 'backfill'):
    # no planB here - we assume it works or we need to know why it fails
    
    # Make a list of files to use from the hycom1 archive.
    # first get a list of all available times
    h_in_dir = Ldir['data'] + 'hycom1/'
    h_list0 = os.listdir(h_in_dir)
    h_list0.sort()
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
    # note that we don't actually write any of the files, because
    # they already exist in LiveOcean_data/hycom1/

if planB == False:
    # process the hycom files
    
    if Ldir['run_type'] == 'forecast':
        # new forecast version
        # copy in the coordinates (assume those from first file work)
        this_h_dict = pickle.load(open(h_in_dir + h_list[0], 'rb'))
        coord_dict = dict()
        for vn in ['lon', 'lat', 'z']:
            coord_dict[vn] = this_h_dict[vn]
        c_out_dir = Ldir['LOogf_fd']
        pickle.dump(coord_dict, open(c_out_dir + 'coord_dict.p', 'wb'))
        
    elif Ldir['run_type'] == 'backfill':
        # backfill
        c_out_dir = Ldir['LOogf_fd']
        coords_dict = pickle.load(open(h_in_dir + 'coords_dict.p', 'rb'))
        coord_dict = coords_dict['91.2']
        pickle.dump(coord_dict, open(c_out_dir + 'coord_dict.p', 'wb'))

    #%% filter in time
    fh_dir = Ldir['LOogf_fd']
    Ofun.time_filter(h_in_dir, h_list, fh_dir, Ldir)

    #%% extrapolate
    lon, lat, z, L, M, N, X, Y = Ofun.get_coords(fh_dir)
    a = os.listdir(fh_dir)
    a.sort()
    aa = [item for item in a if item[:2]=='fh']
    for fn in aa:
        print('-Extrapolating ' + fn)
        in_fn = fh_dir + fn
        V = Ofun.get_extrapolated(in_fn, L, M, N, X, Y, lon, lat, z, Ldir,
            add_CTD=add_CTD, fix_NSoG=fix_NSoG)
        pickle.dump(V, open(fh_dir + 'x' + fn, 'wb'))

    # and interpolate to ROMS format
    # get grid and S info
    G = zrfun.get_basic_info(Ldir['grid'] + 'grid.nc', only_G=True)
    S_info_dict = Lfun.csv_to_dict(Ldir['grid'] + 'S_COORDINATE_INFO.csv')
    S = zrfun.get_S(S_info_dict)
    # get list if files to work on
    a = os.listdir(fh_dir)
    a.sort()
    aa = [item for item in a if item[:3]=='xfh']
    # HyCOM grid info
    lon, lat, z, L, M, N, X, Y = Ofun.get_coords(fh_dir)
    # load a dict of hycom fields
    dt_list = []
    count = 0
    c_dict = dict()
    
    if testing == False: # this is a slow step, so omit for testing
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

elif planB == True:
    print('**** Using planB ****')
    ds_today = Ldir['date_string']
    dt_today = datetime.strptime(ds_today, '%Y.%m.%d')
    dt_yesterday = dt_today - timedelta(days=1)
    ds_yesterday = datetime.strftime(dt_yesterday, format='%Y.%m.%d')
    clm_yesterday = (Ldir['LOog'] + 'f' + ds_yesterday + '/'
        + Ldir['frc'] + '/' + 'ocean_clm.nc')
    clm_today = Ldir['LOogf_f'] + 'ocean_clm.nc'
    shutil.copyfile(clm_yesterday, clm_today)
    ds = nc.Dataset(clm_today, 'a')
    ot = ds['ocean_time'][:]
    ot[-1] += 86400
    for tname in ['ocean', 'salt', 'temp', 'v3d', 'v2d', 'zeta']:
        ds[tname + '_time'][:] = ot
    ds.close()

#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)

if testing == False:
    nc_dir = Ldir['LOogf_f']
    if do_bio and (planB==False):
        G = zrfun.get_basic_info(Ldir['grid'] + 'grid.nc', only_G=True)
        Ofun_bio.add_bio(nc_dir, G, add_CTD=add_CTD)
    Ofun_nc.make_ini_file(nc_dir)
    Ofun_nc.make_bry_file(nc_dir)
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
else:
    result_dict['result'] = 'testing'

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)


