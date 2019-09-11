# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:11:16 2016

@author: PM5

Some extra river functions.

"""

import os
import netCDF4 as nc
import pandas as pd
import river_class
import numpy as np
from datetime import datetime
import Lfun # assume path is provided by calling function

ncformat = 'NETCDF3_64BIT_OFFSET' # NETCDF3_CLASSIC'
    
def dt64_to_dt(dt64):
    # convert numpy datetime64 to datetime
    dt = datetime.utcfromtimestamp(dt64.astype('datetime64[ns]').tolist()/1e9)
    return dt
    
def write_to_nc(out_fn, S, df, qt_df_dict, dt_ind):
    
    # get rid of the old version, if it exists
    try:
        os.remove(out_fn)
    except OSError:
        pass # assume error was because the file did not exist
    foo = nc.Dataset(out_fn, 'w', format=ncformat)
    
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
        csw = S['Cs_w']
        dcsw = np.diff(csw)
        v_var[:, count] = dcsw
        count += 1
        # should end up with velocity constant over depth
    v_var.long_name = 'river runoff mass transport vertical profile'
    v_var.requires = "must sum to 1 over s_rho"

    foo.close()
    

