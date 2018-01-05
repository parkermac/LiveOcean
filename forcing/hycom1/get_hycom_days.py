"""
Extract and save data from a sequence of hycom past days.

At home this takes 70-80 seconds per day.

"""

# setup
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

from importlib import reload
import hfun
reload(hfun)

import netCDF4 as nc
import pickle

# specify the output directory
out_dir = Ldir['data'] + 'hycom1/'
Lfun.make_dir(out_dir, clean=False)

if Ldir['env'] == 'pm_mac': # mac version
    testing = True
elif Ldir['env'] == 'pm_fjord': # fjord version
    testing = False

if testing:
    exnum_list = ['90.9']
else:
    #exnum_list = ['90.9', '91.0', '91.1', '91.2']
    exnum_list = ['91.2']

dt_list_new = []

for exnum in exnum_list:
    # open the Dataset for this exnum
    print('\nWorking on exnum ' + exnum)
    fn = 'http://beta.hycom.org/thredds/dodsC/GLBu0.08/expt_' + exnum  
    
    # get time and space coordinates for this exnum
    ds = nc.Dataset(fn)
    coords = hfun.get_coordinates(ds)    
    dt_list = hfun.get_dt_list(ds)
    ds.close()
    
    # Limiting downloads when testing; counter only increments
    # when we actually save a new, valid file.
    counter = 0 
    if testing:
        maxcount = 1 # the code will get this many good files
    else:
        maxcount = len(dt_list)
    
    for nt in range(len(dt_list)):        
    
        # name the file based on date
        dt = dt_list[nt]
        dts = dt.strftime('%Y.%m.%d')
        print('\n  working on %s' % dts)
        sys.stdout.flush()
        out_name = 'h' + dts + '.p'
        out_fn = out_dir + out_name
        
        # Check to see if file exists.
        # (we could also do this -faster?- by comparing to a directory listing)
        if os.path.exists(out_fn)== True:
            print('  file exists for this day')
            sys.stdout.flush()
            # and in this case we don't replace it           
        else:
            if counter < maxcount:
                # get the data
                ds = nc.Dataset(fn)              
                out_dict = hfun.get_hycom_day(ds, nt, coords)                        
                ds.close()
                # and if it is valid, write it to a file
                if out_dict['result'] == 'success':
                    out_dict['dt'] = dt                
                    print('  writing ' + out_name)
                    sys.stdout.flush()
                    pickle.dump(out_dict, open(out_fn, 'wb'))                    
                    dt_list_new.append(dt)
                    counter+=1                    
            else:
                print('  reached maximum file number')
                sys.stdout.flush()
                break
    