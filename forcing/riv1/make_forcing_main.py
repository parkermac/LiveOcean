"""
This is the main program for making the RIV forcing file.

For development run from the command line in spyder as, e.g.:

cd /Users/PM5/Documents/LiveOcean/forcing/riv1

test historical:
run make_forcing_main.py -g cascadia2 -t frc2 -r backfill -d 2013.01.01

test crossing year boundary:
run make_forcing_main.py -g cascadia2 -t frc2 -r backfill -d 2014.12.31

test times more recent than historical files:
run make_forcing_main.py -g cascadia2 -t frc2 -r backfill -d 2016.01.19

test forecast:
run make_forcing_main.py -g cascadia2 -t frc2 -r forecast

NOTE: all use a different grid than that specified
in the ffun.intro() defaults

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
import river_class
reload(river_class)
import river_functions as rivfun
reload(rivfun)

import pandas as pd
import numpy as np
import netCDF4 as nc


#%% ****************** CASE-SPECIFIC CODE *****************

# set up the time index for the record

from datetime import datetime, timedelta
start_time = datetime.now()

# set first and last times to be at noon
dt0 = datetime.strptime(Ldir['date_string'],'%Y.%m.%d') - timedelta(days=2.5)
dt1 = datetime.strptime(Ldir['date_string'],'%Y.%m.%d') + timedelta(days=4.5)
days = (dt0, dt1)

day_list = []
this_day = dt0
while this_day <= dt1:
    day_list.append(this_day)
    this_day += timedelta(days=1)

# pandas Index objects
dt_ind = pd.Index(day_list)
yd_ind = pd.Index(dt_ind.dayofyear)

# seve some info
Info = dict()
Info['run_type'] = Ldir['run_type']
Info['datestring_start'] = dt0.strftime('%Y.%m.%d')
Info['datestring_end'] = dt1.strftime('%Y.%m.%d')

#%% Load a dataframe with info for rivers to get

ri_fn = Ldir['run'] + 'river_info.csv'
df = pd.read_csv(ri_fn, index_col='rname')

#%% associate rivers with ones that have temperature climatology data

df = rivfun.get_tc_rn(df)

#%% step through all rivers

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

qt_df_dict = dict()

for rn in df.index:

    print(10*'#' + ' ' + rn + ' ' + 10*'#')

    # initialize a qt (flow vs. time) DataFrame for this river
    qt_df = pd.DataFrame(index=dt_ind,
                         columns=['clim','his','usgs','ec','nws','final',
                                  'temperature'])

    rs = df.ix[rn] # a Series with info for this river
    riv = river_class.River(rn, rs)

    # Get climatology (using squeeze=True returns a Series)
    clim = pd.read_csv(Ldir['data'] + 'rivers/Data_clim/' + rn + '.csv',
                    header=None, index_col=0, squeeze=True)
    qt_clim_yd = clim.ix[yd_ind] # clip just the needed values

    # Get T climatology (using squeeze=True returns a Series)
    tc_rn = df.ix[rn, 'tc_rn']
    T_clim = pd.read_csv(Ldir['data'] + 'rivers/Data_T_clim/' + tc_rn + '.csv',
                    header=None, index_col=0, squeeze=True)
    T_clim_yd = T_clim.ix[yd_ind] # clip just the needed values

    # start to populate the qt DataFrame
    qt_df['clim'] = pd.Series(index=dt_ind, data=qt_clim_yd.values)

    # Get historical record (a Series)
    his = pd.read_pickle(Ldir['data'] + 'rivers/Data_historical/'
                + rn + '.p')

    if dt1 <= his.index[-1]:
        # fill with historical data if the timing is right
        qt_df['his'] = his.ix[dt_ind]
        qt_df['final'] = qt_df['his']
        print(' filled from historical')
    else:
        # otherwise try (sequentially) to fill from
        # nws, or usgs, or ec
        if pd.notnull(rs.nws) and Ldir['run_type'] == 'forecast':
            riv.get_nws_data()
            if not riv.qt.empty:
                qt_df['nws'] = riv.qt.ix[dt_ind]
                qt_df['final'] = qt_df['nws']
                print(' filled from nws forecast')
        elif pd.notnull(rs.usgs):
            riv.get_usgs_data(days)
            if not riv.qt.empty:
                qt_df['usgs'] = riv.qt.ix[dt_ind]
                qt_df['final'] = qt_df['usgs']
                print(' filled from usgs')
        elif pd.notnull(rs.ec):
            riv.get_ec_data(days)
            if not riv.qt.empty:
                qt_df['ec'] = riv.qt.ix[dt_ind]
                qt_df['final'] = qt_df['ec']
                print(' filled from ec')

    # check results and fill with extrapolation (ffill) or climatology

    if False: # introduce errors for testing
        qt_df.ix[-3:, 'final'] = np.nan

    if ( pd.isnull(qt_df['final'].values).any() and
            not pd.isnull(qt_df['final'].values).all() ):
        qt_df['final'] = qt_df['final'].ffill(axis=0)
        print(' extended by ffill')

    if pd.isnull(qt_df['final'].values).any():
        qt_df['final'] = qt_df['clim']
        print( 'WARNING: missing values: all filled with climatology')

    if (qt_df['final'].values < 0).any():
        qt_df['final'] = qt_df['clim']
        print( 'WARNING: negative values: all filled with climatology')

    if pd.isnull(qt_df['final'].values).any():
        print( '>>>>>>> flow has missing values!! <<<<<<<<<')

    # Temperature data

    qt_df['temperature'] = pd.Series(index=dt_ind, data=T_clim_yd.values)

    if pd.isnull(qt_df['temperature'].values).any():
        print( '>>>>>>> temp has missing values!! <<<<<<<<<')

    # save in the dict
    qt_df_dict[rn] = qt_df


#%% calculations for vertical distribution

S_info_dict = Lfun.csv_to_dict(Ldir['grids']+Ldir['gridname']+'/S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)

#%% save the output to NetCDF

out_fn = (Ldir['LOogf_f'] + 'rivers.nc')
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
foo = nc.Dataset(out_fn, 'w', format='NETCDF3_CLASSIC')

nriv = len(df)
N = S['N']
ndt = len(dt_ind)
SL = 50 # max string length for river names

foo.createDimension('river', nriv)
foo.createDimension('s_rho', N)
foo.createDimension('river_time', ndt)
foo.createDimension('slen', SL)

v_var = foo.createVariable('river', float, ('river'))
v_var[:] = np.arange(1, nriv+1)
v_var.long_name = 'river runoff identification number'

v_var = foo.createVariable('river_name', 'c', ('slen', 'river'))
rr = 0
for rn in df.index:
    print('- ' + rn)
    cc = 0
    for ch in rn:
        v_var[cc, rr] = ch
        cc += 1
    rr += 1

def dt64_to_dt(dt64):
    # convert numpy datetime64 to datetime
    dt = datetime.utcfromtimestamp(dt64.astype('datetime64[ns]').tolist()/1e9)
    return dt

v_var = foo.createVariable('river_time', float, ('river_time'))
count = 0
for item in dt_ind.values:
    item_dt = dt64_to_dt(item)
    v_var[count] = Lfun.datetime_to_modtime(item_dt)
    count += 1
v_var.long_name = 'river runoff time'
v_var.units = "seconds since 1970-01-01 00:00:00"

v_var = foo.createVariable('river_direction', float, ('river'))
count = 0
for rn in df.index:
    v_var[count] = df.ix[rn, 'idir']
    count += 1
v_var.long_name = 'river runoff direction'
v_var.flag_values = "0, 1"
v_var.flag_meanings = "flow across u-face, flow across v-face"
v_varLwSrc_True = "flag not used"

v_var = foo.createVariable('river_Xposition', float, ('river'))
count = 0
for rn in df.index:
    if df.ix[rn, 'idir'] == 0:
        v_var[count] = df.ix[rn, 'col_py'] + 1
    elif df.ix[rn, 'idir'] == 1:
        v_var[count] = df.ix[rn, 'col_py']
    count += 1
v_var.long_name = 'river XI-position'
v_var.LuvSrc_True_meaning = "i point index of U or V face source/sink"
v_var.LwSrc_True_meaning = "i point index of RHO center source/sink" ;

v_var = foo.createVariable('river_Eposition', float, ('river'))
count = 0
for rn in df.index:
    if df.ix[rn, 'idir'] == 0:
        v_var[count] = df.ix[rn, 'row_py']
    if df.ix[rn, 'idir'] == 1:
        v_var[count] = df.ix[rn, 'row_py'] + 1
    count += 1
v_var.long_name = 'river ETA-position'
v_var.LuvSrc_True_meaning = "j point index of U or V face source/sink"
v_var.LwSrc_True_meaning = "j point index of RHO center source/sink" ;

v_var = foo.createVariable('river_transport', float, ('river_time', 'river'))
count = 0
for rn in df.index:
    qt_df = qt_df_dict[rn]
    flow = qt_df['final'].values
    v_var[:, count] = flow * df.ix[rn, 'isign']
    count += 1
v_var.long_name = 'river runoff vertically integrated mass transport'
v_var.positive = "LuvSrc=T flow in positive u,v direction, LwSrc=T flow into RHO-cell"
v_var.negative = "LuvSrc=T flow in negative u,v direction, LwSrc=T flow out of RHO-cell"
v_var.time = "river_time"
v_var.units = "meter3 second-1"

v_var = foo.createVariable('river_temp', float, ('river_time', 's_rho', 'river'))
count = 0
for rn in df.index:
    qt_df = qt_df_dict[rn]
    for nn in range(N):
        v_var[:, nn, count] = qt_df['temperature'].values
    count += 1
v_var.long_name = 'river runoff potential temperature'
v_var.time = "river_time"
v_var.units = "Celsius"

v_var = foo.createVariable('river_salt', float, ('river_time', 's_rho', 'river'))
count = 0
for rn in df.index:
    for nn in range(N):
        v_var[:, nn, count] = np.zeros(ndt)
    count += 1
v_var.long_name = 'river runoff salinity'
v_var.time = "river_time"
v_var.units = "psu"

v_var = foo.createVariable('river_Vshape', float, ('s_rho', 'river'))
count = 0
for rn in df.index:
    # copied from old matlab code, and simplified
    #   %linear decay from surface value to 0, fixed by sng 7/2011
    v_var[:, count] = np.linspace(0,2/N,N)
    count += 1
v_var.long_name = 'river runoff mass transport vertical profile'
v_var.requires = "must sum to 1 over s_rho"

foo.close()

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
