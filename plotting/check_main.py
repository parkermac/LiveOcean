"""
Check the output of make_forcing_main.py.
"""
# specify which forcing this code is for
which_force = 'ocn'

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
import zfun

Info = Lfun.csv_to_dict(Ldir['LOo'] + 'current_Info/'
    + which_force + '/Info_for_main.csv')
    
# override date
Info['date_string'] = '2015.03.23'
Info['f_string'] = 'f' + Info['date_string']
    
# define the output location
nc_dir = Ldir['LOo'] + Ldir['gtag'] + '/' + Info['f_string'] + '/' + which_force + '/Data/'

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

import netCDF4 as nc

# for reference, here are the choices for fld_name
# fld_list = ['ssh', 's3d', 't3d', 'u3d', 'v3d']
fld_name = 's3d'

fn = nc_dir + fld_name + '.nc'
ds = nc.Dataset(fn)

for vn in ds.variables:
    print(ds.variables[vn])
        
fld = ds.variables[fld_name][:]
flde = ds.variables[fld_name + '_extrap'][:]
fldf = ds.variables[fld_name + '_filt'][:]
lon = ds.variables['lon'][:]
lat = ds.variables['lat'][:]
dt = ds.variables['dt'][:]

# get time series
nr = 80; nc = 40
if fld_name == 'ssh':
    fld_top = fld[0,:,:].squeeze()
    flde_top = flde[0,:,:].squeeze()
    fldf_top = fldf[0,:,:].squeeze()    
    fs = fld[:,nr,nc].squeeze()
    fse = flde[:,nr,nc].squeeze()
    fsf = fldf[:,nr,nc].squeeze()
else:
    nz = 10
    fld_top = fld[0,nz,:,:].squeeze()
    flde_top = flde[0,nz,:,:].squeeze()
    fldf_top = fldf[0,nz,:,:].squeeze() 
    fs = fld[:,nz,nr,nc].squeeze()
    fse = flde[:,nz,nr,nc].squeeze()
    fsf = fldf[:,nz,nr,nc].squeeze()
        
# get the coastline
coast_file = Ldir['data'] + 'coast/pnw_coast_combined.mat'
import scipy.io
cmat = scipy.io.loadmat(coast_file)
clon = cmat['lon'].squeeze()
clat = cmat['lat'].squeeze()

# PLOTTING
import matplotlib.pyplot as plt
plt.close()

fig = plt.figure(figsize=(20, 20))
cmap = plt.get_cmap(name='rainbow')

fl = dict()
fl['fld_top'] = fld_top
fl['flde_top'] = flde_top
fl['fldf_top'] = fldf_top
fll = ['fld_top','flde_top','fldf_top']
for ii in range(3):
    ax0 = fig.add_subplot(2,3,ii+1)
    cs0 = ax0.pcolormesh(lon, lat, fl[fll[ii]], cmap=cmap)# vmin=30, vmax=33)
    ax0.plot(lon[nr,nc], lat[nr,nc],'*r')
    # add coastline
    ax0.plot(clon, clat, '-k', linewidth=.5)
    # limits, scaling, and labels
    ax0.axis([lon.min(), lon.max(), lat.min(), lat.max()])
    zfun.dar(ax0)
    ax0.set_xlabel('Longitude')
    ax0.set_ylabel('Latitude')
    ax0.set_title(fld_name + ' ' + str(dt[0]/86400.))
    fig.colorbar(cs0)

ax1 = fig.add_subplot(212)
ax1.plot((dt - dt[0])/86400, fs,'-*b', linewidth=3 )
ax1.plot((dt - dt[0])/86400, fse,'-og', linewidth=3 )
ax1.plot((dt - dt[0])/86400, fsf,'-sr', linewidth=3 )
plt.show()
