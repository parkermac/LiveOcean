# -*- coding: utf-8 -*-
"""
This is the main program for making the OCN forcing file.

It is designed to work with the new hycom1 archive of extracted files,
as well as regular forecasts.

*******************************

To run from the command line in LiveOcean/driver/:
    
./driver_forcing1.sh -g cas1 -t base -r backfill -0 20130101 -1 20130101

To test in python on mac:

cd /Users/PM5/Documents/LiveOcean/forcing/ocn1

This is a normal case (no gaps)
run make_forcing_main.py -g cas1 -t base -r backfill -d 2013.01.01

This has a huge gap, but still works:
run make_forcing_main.py -g cas1 -t base -r backfill -d 2012.12.01

Or try a forecast
This has a huge gap, but still works:
run make_forcing_main.py -g cas1 -t base -r forecast

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
#import shutil
import pickle

import zrfun
import Ofun
from importlib import reload
reload(Ofun)

#import Ofun_nc

start_time = datetime.now()

# get grid and S info
[G] = zrfun.get_basic_info(Ldir['grid'] + 'grid.nc', getS=False, getT=False)
S_info_dict = Lfun.csv_to_dict(Ldir['grid'] + 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)

vnl_full = ['ssh','s3d','t3d','u3d','v3d']
exnum = '91.2'

if Ldir['run_type'] == 'forecast':    
    h_out_dir = out_dir = Ldir['LOogf_fd']      
    print('** START getting catalog')
    # create a list of url's of the preprocessed HYCOM files for this forecast
    fn_list = Ofun.get_hycom_file_list(exnum)
    # get a selection of the raw list (e.g. one per day)
    varf_dict, dt_list = Ofun.get_varf_dict(fn_list, Ldir)
    var_list = list(varf_dict.keys())    
    vnl_dict = {'ssh':['ssh'], 'ts3z':['s3d','t3d'], 'uv3z':['u3d','v3d']}
    #get the data and pack it in pickle files
    for vns in var_list:
        this_list = varf_dict[vns]
        for fn in this_list:# [:2]: # debugging
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
       
elif Ldir['run_type'] == 'backfill':
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

# copy in the coordinates (assume those from hycom1 work)
exnum1 = '91.2'
c_in_dir = Ldir['data'] + 'hycom1/'
c_out_dir = out_dir = Ldir['LOogf_fd']
coords_dict = pickle.load(open(c_in_dir + 'coords_dict.p', 'rb'))
coord_dict = dict()
for vn in ['lon', 'lat', 'z']:
    coord_dict[vn] = coords_dict[exnum1][vn]
pickle.dump(coord_dict, open(c_out_dir + 'coord_dict.p', 'wb'))    
        
# filter in time
fh_dir = Ldir['LOogf_fd']
Ofun.time_filter(h_in_dir, h_list, fh_dir, Ldir)


#%% extrapolate
lon, lat, z, L, M, N, X, Y = Ofun.get_coords(fh_dir)
a = os.listdir(fh_dir)
aa = [item for item in a if item[:2]=='fh']
for fn in aa:    
    print('\nExtrapolating ' + fn)    
    in_fn = fh_dir + fn
    V = Ofun.get_extrapolated(in_fn, L, M, N, X, Y, z)
    pickle.dump(V, open(fh_dir + 'x' + fn, 'wb'))

# and interpolate to ROMS format
#Ofun.interpolate_to_roms(data_dir, G, S)
#
## write in ROMS format
#Ofun.write_roms_files(data_dir, out_dir)
#
#
#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)

# code for checking temporary results
fh_list0 = os.listdir(fh_dir)
fh_list = [item for item in fh_list0 if item[:2]=='fh']
a = pickle.load(open(fh_dir + fh_list[0], 'rb'))
dt0 = a['dt']
a = pickle.load(open(fh_dir + fh_list[-1], 'rb'))
dt1 = a['dt']

result_dict['var_start_time'] = dt0.strftime(time_format)
result_dict['var_end_time'] = dt1.strftime(time_format)

if os.path.isfile(fh_dir + fh_list[-1]):
   result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)


