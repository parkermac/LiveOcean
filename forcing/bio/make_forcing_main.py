"""
This is the main program for making the BIO forcing file.
"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************
import shutil
fn_list = []
fn_list.append((Ldir['LOo'] + Ldir['gtag'] +
    '/f' + Ldir['date_string'] + '/ocn/ocean_clm'))
fn_list.append((Ldir['LOo'] + Ldir['gtag'] +
    '/f' + Ldir['date_string'] + '/ocn/ocean_bry'))
fn_list.append((Ldir['LOo'] + Ldir['gtag'] +
    '/f' + Ldir['date_string'] + '/ocn/ocean_ini'))
fn_list.append((Ldir['LOo'] + Ldir['gtag'] +
    '/f' + Ldir['date_string'] + '/riv/rivers'))
for fn in fn_list:
    fn_orig = (fn + '.nc')
    fn_for_bio = (fn + '_bio.nc')
    # first get rid of any existing copies of the forcing that
    # might have bio variables, because the worker needs them to be absent
    try:
        os.remove(fn_for_bio)
    except OSError:
        pass
        # assume the file did not exist
    # then make copies with new names
    shutil.copyfile(fn_orig,fn_for_bio) 

# ******************************************************* 

# run the code to create the forcing files
Lfun.run_worker(Ldir)

from datetime import datetime
print('MAIN end time = ' + str(datetime.now()))
