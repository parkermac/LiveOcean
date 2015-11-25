"""
Plots the result of the river forcing code.
"""

# setup
import os
import sys   
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

rpath = Ldir['LO'] + 'forcing/riv/'
if rpath not in sys.path:
    sys.path.append(rpath)
import Rfun
import zfun

# where is the output
indir0 = Ldir['LOo'] + Ldir['gtag'] + '/'

# make a list of date strings
from datetime import datetime, timedelta
ds_list = []
dt0 = datetime(2014, 11, 20)
dt1 = datetime(2015, 1, 27)
dt = dt0
while dt <= dt1:
    ds_list.append(dt.strftime('%Y.%m.%d'))
    dt = dt + timedelta(1)
    
# get the list of river names
# (much simpler version than in LO river code)
fn = Ldir['run'] + 'rname_list.txt'
ff = open(fn, 'r')
rn_list = ff.read().splitlines()
ff.close()
# swap some names
rn_list[rn_list.index('duwamish')] = 'green'
rn_list[rn_list.index('hammahamma')] = 'hamma'

# initialize some information lists
code_list = []
scale_list = []
# step through the river names    
for rn in rn_list:   
    # this function reads a file that has scale factors
    # linked to USGS numbers
    code_ecology, scale_factor_ecology = Rfun.get_river_code_ecology(rn, Ldir)
    code_list.append(code_ecology)
    scale_list.append(scale_factor_ecology)
# make dicts
code_dict = dict(zip(rn_list, code_list))
scale_dict = dict(zip(rn_list, scale_list))

# simpler and more correct version
Rdf, Rnames = Rfun.get_rivers_dataframe(Ldir)

# get the climatology and plot it
import pandas as pd
# for the flow
clim_fn = (Ldir['LO'] + 'forcing/riv/' +
    'river_climatology/Qclim.csv')
clim_df = pd.read_csv(clim_fn, index_col='yearday')
# and for temperature
tclim_fn = (Ldir['LO'] + 'forcing/riv/' +
    'river_climatology/Tclim.csv')
tclim_df = pd.read_csv(tclim_fn, index_col='yearday')
#
# the two data frames have yearday (1:366) as index, and river names for columns

import matplotlib.pyplot as plt
plt.close()
NR = 4
NC = 4
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,8))
count = 0
import numpy as np
for rn in rn_list:
    if rn != 'skagit_south':
        nr = np.divide(count, NR)
        nc = np.mod(count, NC)
        x = clim_df.index.values
        y = clim_df[rn].values * scale_dict[rn]
        ax = axes[nr, nc]
        ax.plot(x, y, '-b')
        count = count + 1
        ax.set_xlim(0, 366)
        ax.text(.1, .83, rn.capitalize(), transform=ax.transAxes)
        ax.text(.1, .7, code_dict[rn], transform=ax.transAxes)
        if nr == NR-1:
            ax.set_xlabel('Yearday')
        
# now add line segments from a forecast
import netCDF4 as nc4
ds_list_short = ds_list #[ds_list[0]]
for ds in ds_list_short:
    riv_nc_fn = (Ldir['LOo'] + Ldir['gtag'] + '/f' + ds +
        '/riv/rivers.nc')
    rds = nc4.Dataset(riv_nc_fn,'r')
    rt_vec = rds.variables['river_time'][:]
    rq = rds.variables['river_transport'][:]
    rdt_list = []   
    for rt in rt_vec:
        rdt_list.append(Lfun.modtime_to_datetime(rt))
    rdf = pd.DataFrame(rq, index=rdt_list, columns=rn_list)
    count = 0
    for rn in rn_list:
        if rn != 'skagit_south':
            nr = np.divide(count, NR)
            nc = np.mod(count, NC)
            x = rdf.index.dayofyear
            y = rdf[rn].values
            ax = axes[nr, nc]
            ax.plot(x, np.abs(y), '.r')
            count = count + 1
    

#zfun.ncd(rds)  
    
    
plt.show()
