"""
This is the main program for reorganizing the files from a daily forecast.
It is designed to work only on a daily forecast folder that has history
files 1-73 in it, and will exit if that is not what it finds.

It creates new folders for the subsequent two days and populates them
with copies of the history files:
day 2: 25-49 become 1-25
day3: 49-73 become 1-25

"""

import os, sys
sys.path.append(os.path.abspath('../'))
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************
from datetime import datetime, timedelta
import shutil

start_time = datetime.now()

ds0 = Ldir['date_string']
dt0 = datetime.strptime(ds0, '%Y.%m.%d')
dt1 = dt0 + timedelta(days=1)
dt2 = dt0 + timedelta(days=2)
ds1 = dt1.strftime('%Y.%m.%d')
ds2 = dt2.strftime('%Y.%m.%d')

print(' - Distributing forecast files for ' + ds0)
f_string0 = 'f' + ds0
f_string1 = 'f' + ds1
f_string2 = 'f' + ds2

dir0 = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string0 + '/'
dir1 = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string1 + '/'
dir2 = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string2 + '/'

fn_list_raw = os.listdir(dir0)
fn_list = [h for h in fn_list_raw if ('ocean_his' in h) and ('.nc' in h)]
fn_list.sort()
NF = len(fn_list)
if NF != 73:
    print('Expecting 73 history files but only %d found.' % (NF))
    sys.exit()
    
testing = False

result = 'success'
try:
    Lfun.make_dir(dir1, clean=True)
    Lfun.make_dir(dir2, clean=True)
    for ii in range(1,26):
        h0 = ('0000' + str(ii+24))[-4:]
        h1 = ('0000' + str(ii))[-4:]
        fn0 = dir0 + 'ocean_his_' + h0 + '.nc'
        fn1 = dir1 + 'ocean_his_' + h1 + '.nc'
        if testing:
            print('\n'+fn0)
            print(fn1)
        else:
            shutil.copyfile(fn0, fn1)
    for ii in range(1,26):
        h0 = ('0000' + str(ii+48))[-4:]
        h2 = ('0000' + str(ii))[-4:]
        fn0 = dir0 + 'ocean_his_' + h0 + '.nc'
        fn2 = dir2 + 'ocean_his_' + h2 + '.nc'
        if testing:
            print('\n'+fn0)
            print(fn2)
        else:
            shutil.copyfile(fn0, fn2)
except Exception as e:
    result = 'fail'
    print(e)

# ======================================================================

#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)
result_dict['result'] = result

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)
    

