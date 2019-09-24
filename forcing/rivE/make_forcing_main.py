"""
This is the main program for making the RIV forcing file.

test command:

run make_forcing_main.py -g aestus3 -t v1 -r backfill -d 2013.01.01

NOTE: this is designed for hand-manipulation of the river forcing,
such as making an artifical exchange flow for the aestus3 grid
for Elizabeth.

"""

import os; import sys
sys.path.append(os.path.abspath('../'))

import forcing_functions as ffun
Ldir, Lfun = ffun.intro()
import zrfun

import pandas as pd
import netCDF4 as nc
import numpy as np

#%% ****************** CASE-SPECIFIC CODE *****************

ncformat = 'NETCDF3_64BIT_OFFSET' # NETCDF3_CLASSIC'

# this code is so specific to a single grid that I will enfore the name here
if Ldir['gridname'] != 'aestus3':
    sys.exit()

testing = True
if testing:
    import matplotlib.pyplot as plt
    plt.close('all')
    
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
# save some info
Info = dict()
Info['run_type'] = Ldir['run_type']
Info['datestring_start'] = dt0.strftime(dsf)
Info['datestring_end'] = dt1.strftime(dsf)

#%% get dict S
S_info_dict = Lfun.csv_to_dict(Ldir['grid'] + 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)

#%% name output file
out_fn = (Ldir['LOogf_f'] + 'rivers.nc')

# load grid and find info for source locations
dsg = nc.Dataset(Ldir['grid'] + 'grid.nc')

# find the first rho index from the East that has unmasked cells
mask_rho = dsg['mask_rho'][:]
mask_u = dsg['mask_u'][:]
h = dsg['h'][:]
NR_rho, NC_rho = mask_rho.shape
all_masked = True
cc = NC_rho-1
while all_masked and (cc >= 0):
    this_col = mask_rho[:,cc]
    if (this_col == 1).any():
        print(cc)
        break
    else:
        cc -= 1
this_h = h[:,cc]
this_lat = dsg['lat_rho'][:,cc]
this_dy = 1/dsg['pn'][:,cc]
mask = this_col == 1
hh = this_h[mask].data
yy = this_lat[mask].data
YY = yy * np.ones((S['N'],1))
dy = this_dy[mask].data
zz_rho, zz_w = zrfun.get_z(hh, 0*hh, S)
dz = np.diff(zz_w, axis=0)
da = dz * dy # area of each grid box (on the rho grid)

# set transports and their spatial distribution
H1 = 5 # depth of no motion (m)
h1 = np.minimum(H1,hh)
parabola_1 = -(zz_rho + h1) * (-zz_rho + h1)
parabola_2 = -(zz_rho + hh) * (zz_rho + h1)

# set target transports
Q1 = -15000
Q2 = 14000

h1_mask = zz_rho > -h1

# sum of weights
sp_1 = np.sum(parabola_1[h1_mask])
sp_2 = np.sum(parabola_2[~h1_mask])

q1 = Q1 * parabola_1 / sp_1
q2 = Q2 * parabola_2 / sp_2

q = q1.copy()
q[~h1_mask] = q2[~h1_mask]
print(q.sum())

u1 = q1/da
u2 = q2/da

u = q/da

if testing:
    fig = plt.figure(figsize=(14,8))
    
    ax = fig.add_subplot(221)
    ax.plot(yy, -hh, 'og')
    ax.set_ylim(-hh.max()-1, 1)
    ax.grid(True)
    ax.plot(YY[h1_mask], zz_rho[h1_mask], '.b')
    ax.plot(YY[~h1_mask], zz_rho[~h1_mask], '.r')
    ax.set_title('z_rho of two layers')
    
    ax = fig.add_subplot(222)
    ax.plot(parabola_1[h1_mask], zz_rho[h1_mask], '.b')
    ax.plot(parabola_2[~h1_mask], zz_rho[~h1_mask], '.r')
    ax.set_title('parabola shapes')
    
    ax = fig.add_subplot(223)
    ax.plot(q1[h1_mask], zz_rho[h1_mask], '.b')
    ax.plot(q2[~h1_mask], zz_rho[~h1_mask], '.r')
    ax.set_title('transports (m3/s)')
    
    ax = fig.add_subplot(224)
    ax.plot(u1[h1_mask], zz_rho[h1_mask], 'ob')
    ax.plot(u2[~h1_mask], zz_rho[~h1_mask], 'or')
    ax.plot(u, zz_rho,'.c')
    ax.set_title('velocity (m/s)')
    
    plt.show()
    
# write to NetCDF output

# recombine things so that the full boundary condition is
# remapped into NR "rivers" that span the estuary mouth
nriv = len(hh)
ndt = len(day_list)
river_direction = 0 * np.ones(nriv)
col_py = np.ones(nriv, dtype=int) * cc
row_py = np.argwhere(mask)

# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
foo = nc.Dataset(out_fn, 'w', format=ncformat)

N = S['N']

foo.createDimension('river', nriv)
foo.createDimension('s_rho', N)
foo.createDimension('river_time', ndt)

v_var = foo.createVariable('river', float, ('river'))
v_var[:] = np.arange(1, nriv+1)
v_var.long_name = 'river runoff identification number'

v_var = foo.createVariable('river_time', float, ('river_time'))
count = 0
for item in day_list:
    v_var[count] = Lfun.datetime_to_modtime(item)
    count += 1
v_var.long_name = 'river runoff time'
v_var.units = "seconds since 1970-01-01 00:00:00"

v_var = foo.createVariable('river_direction', float, ('river'))
for rr in range(nriv):
    v_var[rr] = river_direction[rr]
v_var.long_name = 'river runoff direction'
v_var.flag_values = "0, 1"
v_var.flag_meanings = "flow across u-face, flow across v-face"
v_varLwSrc_True = "flag not used"

v_var = foo.createVariable('river_Xposition', float, ('river'))
count = 0
for rr in range(nriv):
    if river_direction[rr] == 0:
        v_var[count] = col_py[rr] + 1
    elif river_direction[rr] == 1:
        v_var[count] = col_py[rr]
v_var.long_name = 'river XI-position'
v_var.LuvSrc_True_meaning = "i point index of U or V face source/sink"
v_var.LwSrc_True_meaning = "i point index of RHO center source/sink" ;

v_var = foo.createVariable('river_Eposition', float, ('river'))
for rr in range(nriv):
    if river_direction[rr] == 0:
        v_var[rr] = row_py[rr]
    if river_direction[rr] == 1:
        v_var[rr] = row_py[rr] + 1
v_var.long_name = 'river ETA-position'
v_var.LuvSrc_True_meaning = "j point index of U or V face source/sink"
v_var.LwSrc_True_meaning = "j point index of RHO center source/sink" ;

v_var = foo.createVariable('river_transport', float, ('river_time', 'river'))
for rr in range(nriv):
    this_qr = q[:,rr].sum() * np.ones(ndt)
    v_var[:, count] = this_qr
v_var.long_name = 'river runoff vertically integrated mass transport'
v_var.positive = "LuvSrc=T flow in positive u,v direction, LwSrc=T flow into RHO-cell"
v_var.negative = "LuvSrc=T flow in negative u,v direction, LwSrc=T flow out of RHO-cell"
v_var.time = "river_time"
v_var.units = "meter3 second-1"

v_var = foo.createVariable('river_temp', float, ('river_time', 's_rho', 'river'))
for rr in range(nriv):
    for nn in range(N):
        v_var[:, nn, rr] = 10 * np.ones(ndt)
v_var.long_name = 'river runoff potential temperature'
v_var.time = "river_time"
v_var.units = "Celsius"

v_var = foo.createVariable('river_salt', float, ('river_time', 's_rho', 'river'))
for rr in range(nriv):
    for nn in range(N):
        v_var[:, nn, rr] = np.zeros(ndt)
v_var.long_name = 'river runoff salinity'
v_var.time = "river_time"
v_var.units = "psu"

v_var = foo.createVariable('river_Vshape', float, ('s_rho', 'river'))

for rr in range(nriv):
    v_var[:, rr] = q[:, rr]/q[:,rr].sum()
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
