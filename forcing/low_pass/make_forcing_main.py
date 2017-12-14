"""
This is the main program for making the LOW PASS forcing file.

Performance: about 5 minutes per day on fjord to do
    a low pass of cascadia1_base_lobio1

For testing on my mac run in ipython as

run make_forcing_main.py -d 2017.05.18

run make_forcing_main.py -x lobio1 -d 2017.08.06
"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************
import zfun
import zrfun

from datetime import datetime, timedelta
start_time = datetime.now()

# make input list (full paths)
flist = []
# create the list of history files
if Ldir['run_type'] == 'backfill':
    date_string = Ldir['date_string']
    dt_now = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
    dt_tomorrow = dt_now + timedelta(1)
    dt_yesterday = dt_now - timedelta(1)
    dt_list = [dt_yesterday, dt_now, dt_tomorrow]
    for dt in dt_list:
        date_string = dt.strftime(format='%Y.%m.%d')
        indir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
            '/f' + date_string + '/')
        for ii in range(2,26): # use range(2,26) to use Godin 71 hour filter
            hnum = ('0000' + str(ii))[-4:]
            flist.append(indir + 'ocean_his_' + hnum + '.nc')
    # remove the last item
    flist.pop() # cute
elif Ldir['run_type'] == 'forecast':
    # use the middle day of the last forecast (= yesterday)
    # and today and tomorrow from today's forecast
    date_string = Ldir['date_string']
    dt_now = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
    dt_yesterday = dt_now - timedelta(1)
    dt_list = [dt_yesterday, dt_now]
    for dt in dt_list:
        if dt == dt_yesterday:
            date_string = dt.strftime(format='%Y.%m.%d')
            indir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
                '/f' + date_string + '/')
            for ii in range(2,26): # use range(2,26) to use Godin 71 hour filter
                hnum = ('0000' + str(ii))[-4:]
                flist.append(indir + 'ocean_his_' + hnum + '.nc')
        elif dt == dt_now:
            date_string = dt.strftime(format='%Y.%m.%d')
            indir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
                '/f' + date_string + '/')
            for ii in range(2,49): # use range(2,49) to use Godin 71 hour filter
                hnum = ('0000' + str(ii))[-4:]
                flist.append(indir + 'ocean_his_' + hnum + '.nc')

# make output name (full path)
out_fn = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
    '/f' + Ldir['date_string'] + '/low_passed.nc')
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist

# create the filter
nf = len(flist)
if nf == 71:
    print(' - Using Godin filter')
    filt0 = zfun.godin_shape()
else:
    print(' - Using Hanning filter for list length = ' + str(nf))
    filt0 = zfun.hanning_shape(nf)

# RUN THE FUNCTION
exclude = ['PH', 'ARAG']
zrfun.roms_low_pass(flist, out_fn, filt0, exclude=exclude)

#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)
if os.path.isfile(out_fn):
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)

