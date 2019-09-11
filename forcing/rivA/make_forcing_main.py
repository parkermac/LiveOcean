"""
This is the main program for making the RIV forcing file.

test command:

run make_forcing_main.py -g aestus3 -t v1 -r backfill -d 2013.01.01

NOTE: this is designed for hand-manipulation of the river forcing,
such as making an artifical exchange flow for the aestus3 grid
for Elizabeth.

"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)

from importlib import reload
import forcing_functions as ffun
reload(ffun)
Ldir, Lfun = ffun.intro()
import zrfun
import river_functions as rivfun
reload(rivfun)

import pandas as pd

#%% ****************** CASE-SPECIFIC CODE *****************

# set up the time index for the record
from datetime import datetime, timedelta
start_time = datetime.now()
dsf = '%Y.%m.%d'
# set first and last times to be at noon
dt0 = datetime.strptime(Ldir['date_string'],dsf) - timedelta(days=2.5)
dt1 = datetime.strptime(Ldir['date_string'],dsf) + timedelta(days=4.5)
days = (dt0, dt1)
day_list = []
this_day = dt0
while this_day <= dt1:
    day_list.append(this_day)
    this_day += timedelta(days=1)
# pandas Index objects
dt_ind = pd.Index(day_list)
yd_ind = pd.Index(dt_ind.dayofyear)
# save some info
Info = dict()
Info['run_type'] = Ldir['run_type']
Info['datestring_start'] = dt0.strftime(dsf)
Info['datestring_end'] = dt1.strftime(dsf)

#%% Load a dataframe with info for rivers to get
ri_fn = Ldir['grid'] + 'river_info.csv'
df = pd.read_csv(ri_fn, index_col='rname')

#%% associate rivers with ones that have temperature climatology data
df = rivfun.get_tc_rn(df)

# get the flow and temperature data for these days
qt_df_dict = rivfun.get_qt(df, dt_ind, yd_ind, Ldir, dt1, days)

#%% get dict S
S_info_dict = Lfun.csv_to_dict(Ldir['grid'] + 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)

#%% save the output to NetCDF
out_fn = (Ldir['LOogf_f'] + 'rivers.nc')
rivfun.write_to_nc(out_fn, S, df, qt_df_dict, dt_ind)

#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)
result_dict['var_start_time'] = dt0.strftime(time_format)
result_dict['var_end_time'] = dt1.strftime(time_format)
if os.path.isfile(out_fn):
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)
