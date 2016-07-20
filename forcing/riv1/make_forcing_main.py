"""
This is the main program for making the RIV forcing file.

For development run from the command line in spyder as, e.g.:

cd /Users/PM5/Documents/LiveOcean/forcing/riv1

run make_forcing_main.py -g test -r backfill -d 2015.09.19
run make_forcing_main.py -g test -r backfill -d 2014.12.31
run make_forcing_main.py -g test -r backfill -d 2016.01.19
run make_forcing_main.py -g test -r forecast

to use a different grid than that specified in the ffun.intro() defaults.
"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

from importlib import reload

import zfun
reload(zfun)

import river_class
reload(river_class)

import pandas as pd
#import numpy as np
#import netCDF4 as nc

#%% ****************** CASE-SPECIFIC CODE *****************

# set up the time index for the record

from datetime import datetime, timedelta

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

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

ri_fn = Ldir['run'] + 'river_info.csv'
df = pd.read_csv(ri_fn, index_col='rname')

qt_df_dict = dict()

for rn in ['skagit']:#df.index:

    print(10*'#' + ' ' + rn + ' ' + 10*'#')

    qt_df = pd.DataFrame(index=dt_ind,
                         columns=['clim','his','usgs','ec','nws','final'])

    rs = df.ix[rn] # a series with info for this river
    riv = river_class.River(rn, rs)

    # Get climatology
    # using squeeze=True returns a Series
    clim = pd.read_csv(Ldir['data'] + 'rivers/Data_clim/' + rn + '.csv',
                    header=None, index_col=0, squeeze=True)
    qt_clim_yd = clim.ix[yd_ind]
    qt_df['clim'] = pd.Series(index=dt_ind, data=qt_clim_yd.values)

    # Get historical record (a Series)
    his = pd.read_pickle(Ldir['data'] + 'rivers/Data_historical/'
                + rn + '.p')

    if dt1 <= his.index[-1]:
        # this is a Series
        qt_df['his'] = his.ix[dt_ind]
        qt_df['final'] = qt_df['his']
    else:
        # only one of these will be filled
        if pd.notnull(rs.usgs):
            riv.get_usgs_data(days)
            if not riv.qt.empty:
                qt_df['usgs'] = riv.qt.ix[dt_ind]
                if (pd.isnull(qt_df['usgs'].values).any() and
                    not pd.isnull(qt_df['usgs'].values).all() ):
                    qt_df['usgs'].ffill(axis=0)
        elif pd.notnull(rs.ec):
            riv.get_ec_data(days)
            if not riv.qt.empty:
                qt_df['ec'] = riv.qt.ix[dt_ind]

        if pd.notnull(rs.nws):
            riv.get_nws_data()
            if not riv.qt.empty:
                qt_df['nws'] = riv.qt.ix[dt_ind]
                print('hi')
                qt_df['final'] = qt_df['nws']

    qt_df_dict[rn] = qt_df



#%% calculations for vertical distribution

S_info_dict = Lfun.csv_to_dict(Ldir['run']+'S_COORDINATE_INFO.csv')
S = zfun.get_S(S_info_dict)

#%% write to netcdf

#dimensions:
#	s_rho = 30 ;
#	river = 4  ;
#	river_time = UNLIMITED ; // (0 currently)
#variables:
#	double river(river) ;
#		river:long_name = "river runoff identification number" ;
#	double river_time(river_time) ;
#		river_time:long_name = "river runoff time" ;
#		river_time:units = "days since 2001-01-01 00:00:00" ;
#	double river_direction(river) ;
#		river_direction:long_name = "river runoff direction" ;
#                river_direction:flag_values = "0, 1" ;
#                river_direction:flag_meanings = "flow across u-face, flow across v-face" ;
#                river_direction:LwSrc_True = "flag not used" ;
#	double river_Xposition(river) ;
#		river_Xposition:long_name = "river XI-position" ;
#                river_Xposition:LuvSrc_meaning = "i point index of U or V face source/sink" ;
#                river_Xposition:LwSrc_meaning = "i point index of RHO center source/sink" ;
#	double river_Eposition(river) ;
#		river_Eposition:long_name = "river ETA-position" ;
#                river_Xposition:LuvSrc_True_meaning = "j point index of U or V face source/sink" ;
#                river_Xposition:LwSrc_True_meaning = "j point index of RHO center source/sink" ;
#	double river_transport(river_time, river) ;
#		river_transport:long_name = "river runoff vertically integrated mass transport" ;
#		river_transport:units = "meter3 second-1" ;
#                river_transport:positive = "LuvSrc=T flow in positive u,v direction, LwSrc=T flow into RHO-cell" ;
#                river_transport:negative = "LuvSrc=T flow in negative u,v direction, LwSrc=T flow out of RHO-cell" ;
#		river_transport:time = "river_time" ;
#	double river_Vshape(s_rho, river) ;
#		river_Vshape:long_name = "river runoff mass transport vertical profile" ;
#		river_Vshape:requires = "must sum to 1 over s_rho" ;
#	double river_temp(river_time, s_rho, river) ;
#		river_temp:long_name = "river runoff potential temperature" ;
#		river_temp:units = "Celsius" ;
#		river_temp:time = "river_time" ;
#	double river_salt(river_time, s_rho, river) ;
#		river_salt:long_name = "river runoff salinity" ;
#		river_salt:time = "river_time" ;
#// global attributes:
#		:rivers = "(1) Connecticut River at Hartford, CT, (2) Hudson River at Green Island NY, (3) Penobscot River at Eddington, ME, (4) Delaware River at Trenton NJ " ;
#}

print('MAIN end time = ' + str(datetime.now()))
