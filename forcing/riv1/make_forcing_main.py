"""
This is the main program for making the RIV forcing file.
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

from datetime import datetime, timedelta

Info = dict()
#Info['run_type'] = Ldir['run_type']

Info['run_type'] = 'test'

if Info['run_type'] == 'backfill':
    dt1 = datetime.strptime(Ldir['date_string'],'%Y.%m.%d') + timedelta(3)
    dt0 = datetime.strptime(Ldir['date_string'],'%Y.%m.%d') - timedelta(3)
    days = (dt0, dt1)
    Info['datestring_start'] = dt0.strftime('%Y.%m.%d')
    Info['datestring_end'] = dt1.strftime('%Y.%m.%d')
elif Info['run_type'] == 'forecast':
    days = ()
elif Info['run_type'] == 'test':
    dt0 = datetime(2015,1,1)
    dt1 = datetime(2015,1,30)
    days = (dt0, dt1)
    Info['datestring_start'] = dt0.strftime('%Y.%m.%d')
    Info['datestring_end'] = dt1.strftime('%Y.%m.%d')

#%% Load a dataframe with info for rivers to get

df = pd.read_csv(Ldir['run']+'river_info.csv', index_col='rname')

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
