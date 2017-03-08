"""
Prints information about what hycom archived files are available
"""

# setup
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
import zfun
import netCDF4 as nc
from datetime import datetime, timedelta
exnum_list = ['90.9', '91.0', '91.1', '91.2']
#exnum_list = ['91.1']

for exnum in exnum_list:
    print('\nWorking on ' + exnum)
    fn = 'http://beta.hycom.org/thredds/dodsC/GLBu0.08/expt_' + exnum  
    ds = nc.Dataset(fn)
    if False: # set to True to see what is in the files
        zfun.ncd(ds)
    # get the time in a meaningful format
    t_hycom = ds.variables['time'][:].squeeze()
    tu = ds.variables['time'].units
    print(' time units = ' + tu)  # should be 'hours since 2000-01-01 00:00:00'
    t_origin = ds.variables['time'].time_origin
    dt00 = datetime.strptime(t_origin, '%Y-%m-%d %H:%M:%S')
    dt_list = [] # initialize a list
    for tt in t_hycom:    
        dt_list.append(dt00 + timedelta(tt/24.))
    print(' dt start = ' + datetime.strftime(dt_list[0], '%Y-%m-%d %H:%M:%S'))
    print(' dt end   = ' + datetime.strftime(dt_list[-1], '%Y-%m-%d %H:%M:%S'))
    
    # dt0 = datetime(2014,4,7)
    # dt1 = datetime(2015,3,27)
    #
    # # get full coordinates (vectors for the plaid grid)
    # lon = ds.variables['lon'][:]
    # lat = ds.variables['lat'][:]
    #
    # # find indices of a sub region
    # import numpy as np
    # lon = lon - 360 # convert to -360 to 0 format
    # i0 = np.where(lon > -129)[0][0]
    # i1 = np.where(lon < -121)[-1][-1]
    # j0 = np.where(lat > 41)[0][0]
    # j1 = np.where(lat < 51)[-1][-1]
    #
    # try:
    #     nt0 = dt_list.index(dt0)
    #     print(dt0)
    #     print('nt0 = ' + str(nt0))
    # except:
    #     print('error with dt0!')
    #     nt0 = -1
    # try:
    #     nt1 = dt_list.index(dt1)
    #     print(dt1)
    #     print('nt1 = ' + str(nt1))
    # except:
    #     print('error with dt1!')
    #     nt1 = -1
    #
    # # get some time series
    # v2 = ds.variables['surf_el'][nt0:nt1+1, j0, i0].squeeze()
    # v3 = ds.variables['water_temp'][nt0:nt1+1, 0, j0, i0].squeeze()
    
    ds.close()
