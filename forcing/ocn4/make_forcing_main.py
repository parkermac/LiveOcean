# -*- coding: utf-8 -*-
"""
This is the main program for making the OCN forcing file.

It is designed to work with the new hycom2 archive of extracted files,
as well as regular forecasts.

Performance: on my mac a backfill day for cas6 takes 1 minute.

2017.12.10 I added a planB flag to the forecast case
which copies the clm file from the previous day and
makes the final time a day later.

2018.05.19 I added the add_CTD flag (and the Ofun_CTD module) to add CTD data
on a specified day only.

2019.04.24 I Removed the fix_NSoG code entirely because that experiment did
not produce the desired results.

2019.05.09 Changed the day for which CTD ICs are added to 2016.12.15, but in
Ofun_CTD.get_casts() it is hardwired to look for casts on or after January 2017.

2019.05.20 Added Ofun.get_interpolated_alt() which sped up the program by a factor of 10.

2020.08.13 Added a method to get the HYCOM files using the nco operator "ncks",
which is new essentially Plan A.  I deprecated the "FMRC_best" method to Plan B
because it was failiing about tome time out of four.  Then moved the persistence
backup to Plan C.

*******************************

To run from the command line in LiveOcean/driver/:
    
./driver_forcing2.sh -g cas6 -t v1 -f ocn4 -r backfill -0 20170101 -1 20170101

To test in python on mac:

# standard backfill
run make_forcing_main.py -g cas6 -t v3 -r backfill -d 2017.04.20

# backfill with Salish and coastal estuary IC's from CTD and other info
run make_forcing_main.py -g cas6 -t v1 -r backfill -d 2016.12.15
- the switch to do this is hardwired to a day: 2016.12.15

# today's forecast
run make_forcing_main.py -g cas6 -t v3 -r forecast

"""

import os, sys
sys.path.append(os.path.abspath('../'))
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

#%% ****************** CASE-SPECIFIC CODE *****************

from datetime import datetime, timedelta
import shutil
import pickle
import netCDF4 as nc
import numpy as np
import time

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
planB = False
planC = False
add_CTD = False
do_bio = True

verbose = False

# *** automate when to set add_CTD to True ***
this_dt = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
if this_dt == datetime(2016,12,15):
    print('WARNING: adding CTD data to extrapolation!!')
    add_CTD = True
    
if (Ldir['run_type'] == 'forecast'):
    # this either gets new hycom files, or sets planB to True,
    # and planB may in turn set planC to True
    h_out_dir = Ldir['LOogf_fd']
    Lfun.make_dir(h_out_dir, clean=True)
    # form list of days to get
    nd_f = np.ceil(Ldir['forecast_days'])
    dt0 = this_dt - timedelta(days=2)
    dt1 = this_dt + timedelta(days=int(nd_f) + 2)
    dt_list_full = []
    dtff = dt0
    while dtff <= dt1:
        dt_list_full.append(dtff)
        dtff = dtff + timedelta(days=1)
        
    try:
        # Plan A: use ncks
        # (eventually move to Ofun)
        import subprocess
        def get_hdt_list(fn):
            # function to get the times of a HYCOM ncks extraction as a list of datetimes
            ds = nc.Dataset(fn)
            # get info for the forecast
            t = ds['time'][:]
            if isinstance(t, np.ma.MaskedArray):
                th_vec = t.data
            else:
                th_vec = t
            tu = ds['time'].units
            # e.g. 'hours since 2018-11-20 12:00:00.000 UTC'
            ymd = tu.split()[2]
            hmss = tu.split()[3]
            hms = hmss.split('.')[0]
            hycom_dt0 = datetime.strptime(ymd + ' ' + hms, '%Y-%m-%d %H:%M:%S')
            hdt_list = []
            for th in th_vec:
                this_dt = hycom_dt0 + timedelta(days=(th/24))
                hdt_list.append(this_dt)
            ds.close()
            return hdt_list
        # specify spatial limits
        sys.path.append(os.path.abspath('../hycom2/'))
        import hfun
        aa = hfun.aa
        north = aa[3]
        south = aa[2]
        west = aa[0] + 360
        east = aa[1] + 360
        # name full output file
        full_fn_out = h_out_dir + 'forecast_ncks.nc'
        # time limits
        dstr0 = dt0.strftime('%Y-%m-%dT00:00') 
        dstr1 = dt1.strftime('%Y-%m-%dT00:00')
        # use subprocess.call() to execute the ncks command
        vstr = 'surf_el,water_temp,salinity,water_u,water_v,depth'
        cmd_list = ['ncks',
            '-d', 'time,'+dstr0+','+dstr1+',8',
            '-d', 'lon,'+str(west)+'.,'+str(east)+'.,1',
            '-d','lat,'+str(south)+'.,'+str(north)+'.,1',
            '-v',vstr, 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/latest',
            '-4', '-O', full_fn_out]
        # run ncks
        tt0 = time.time()
        ret = subprocess.call(cmd_list)
        print('Time to get full file using ncks = %0.2f sec' % (time.time()-tt0))
        print('Return code = ' + str(ret) + ' (0=success)')
        # check output
        # - times
        hdt_list = get_hdt_list(full_fn_out)
        if verbose:
            print('\nTarget time range = ' + dstr0 + ' to ' + dstr1)
            for hdt in hdt_list:
                print(' - Actual time = ' + hdt.strftime('%Y-%m-%d-T00:00:00Z'))
            # next split it into individual files as expected by the later processing
            print('')
            print('Split up into separate output files:')
        NT = len(hdt_list)
        for ii in range(NT):
            hdt = hdt_list[ii]
            iis = str(ii)
            fn_out = h_out_dir + 'h'+ hdt.strftime('%Y.%m.%d') + '.nc'
            cmd_list = ['ncks',
                '-d', 'time,'+iis+','+iis,
                '-O', full_fn_out, fn_out]
            ret = subprocess.call(cmd_list)
            this_hdt = get_hdt_list(fn_out)[0]
            if verbose:
                print(fn_out + ': actual time = ' + str(this_hdt))
    except:
        print('- error getting forecast files using ncks')
        planB = True
    
    if planB == True:
        # Plan B: use fmrc one day at a time
        try:
            # get HYCOM files for each day, one day at a time using fmrc
            for dtff in dt_list_full:
                data_fn_out =  h_out_dir + 'h' + dtff.strftime('%Y.%m.%d')+ '.nc'
                if verbose:
                    print('\n' + data_fn_out)
                sys.stdout.flush()
                # get hycom forecast data from the web, and save it in the file "data_fn_out".
                Ofun.get_data_oneday(dtff, data_fn_out)
        except:
            print('- error getting forecast files using fmrc')
            planC = True
       
elif (Ldir['run_type'] == 'backfill'):
    # Make a list of files to use from the hycom2 archive.
    # no planB here - we assume it works or we need to know why it fails    
    h_out_dir = Ldir['LOogf_fd']
    Lfun.make_dir(h_out_dir, clean=True)
    # get the list of raw extractions
    hnc_short_list = Ofun.get_hnc_short_list(this_dt, Ldir)
    # step through those days and convert them to the same format
    # of pickled dicts as used by the forecast
    for fn in hnc_short_list:
        a = Ofun.convert_extraction_oneday(fn) # UNTESTED 2020.04.22
        # a = Ofun.convert_extraction(fn, 0) # function deleted
        dts = datetime.strftime(a['dt'], '%Y.%m.%d')
        out_fn = h_out_dir + 'h' + dts + '.p'
        pickle.dump(a, open(out_fn, 'wb'))
    h_in_dir = h_out_dir
    h_list0 = os.listdir(h_in_dir)
    h_list0.sort()
    h_list = [item for item in h_list0 if item[0] == 'h']
    
if planC == False:
    # process the hycom files
    
    # convert the format
    for dtff in dt_list_full:
        data_fn_out =  h_out_dir + 'h' + dtff.strftime('%Y.%m.%d')+ '.nc'
        a = Ofun.convert_extraction_oneday(data_fn_out)
        out_fn = h_out_dir + 'h' + dtff.strftime('%Y.%m.%d') + '.p'
        pickle.dump(a, open(out_fn, 'wb'))
        sys.stdout.flush()
    h_in_dir = h_out_dir
    h_list0 = os.listdir(h_in_dir)
    h_list0.sort()
    h_list = [item for item in h_list0 if item[0] == 'h' and item[-2:] == '.p']
    
    # copy in the coordinates (assume those from first file work)
    this_h_dict = pickle.load(open(h_in_dir + h_list[0], 'rb'))
    coord_dict = dict()
    for vn in ['lon', 'lat', 'z']:
        coord_dict[vn] = this_h_dict[vn]
    c_out_dir = Ldir['LOogf_fd']
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
            add_CTD=add_CTD)
        pickle.dump(V, open(fh_dir + 'x' + fn, 'wb'))

    # and interpolate to ROMS format
    # get grid and S info
    G = zrfun.get_basic_info(Ldir['grid'] + 'grid.nc', only_G=True)
    S_info_dict = Lfun.csv_to_dict(Ldir['grid'] + 'S_COORDINATE_INFO.csv')
    S = zrfun.get_S(S_info_dict)
    # get list of files to work on
    a = os.listdir(fh_dir)
    a.sort()
    aa = [item for item in a if item[:3]=='xfh']
    # HyCOM grid info
    lon, lat, z, L, M, N, X, Y = Ofun.get_coords(fh_dir)
    # load a dict of hycom fields
    dt_list = []
    count = 0
    c_dict = dict()

    zinds = Ofun.get_zinds(G['h'], S, z)
    for fn in aa:
        print('-Interpolating ' + fn + ' to ROMS grid')
        in_fn = fh_dir + fn
        b = pickle.load(open(in_fn, 'rb'))
        dt_list.append(b['dt'])
        c = Ofun.get_interpolated_alt(G, S, b, lon, lat, z, N, zinds)
        c_dict[count] = c
        count += 1
    #%% Write to ROMS forcing files
    nc_dir = Ldir['LOogf_f']
    tt0 = time.time()
    Ofun_nc.make_clm_file(Ldir, nc_dir, fh_dir, c_dict, dt_list, S, G)
    print(' --write ocean_clm.nc took %0.1f seconds' % (time.time() - tt0))
    
elif planC == True:
    print('**** Using planC ****')
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

nc_dir = Ldir['LOogf_f']
if do_bio and (planC==False):
    tt0 = time.time()
    G = zrfun.get_basic_info(Ldir['grid'] + 'grid.nc', only_G=True)
    Ofun_bio.add_bio(nc_dir, G, add_CTD=add_CTD)
    print(' --add bio took %0.1f seconds' % (time.time() - tt0))
    
tt0 = time.time()
Ofun_nc.make_ini_file(nc_dir)
Ofun_nc.make_bry_file(nc_dir)
print(' --write ocean_ini and _bry.nc took %0.1f seconds' % (time.time() - tt0))
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


