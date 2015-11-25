"""
Check the boundary and clim forcing files files.
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

Info = dict()
# override date
Info['date_string'] = '2015.03.23'
Info['f_string'] = 'f' + Info['date_string']
    

f_dir = (Ldir['LOo'] + Ldir['gtag'] + '/' + Info['f_string'] + '/'
    + which_force + '/')
       
r_dir = (Ldir['roms'] + 'output/' + Ldir['gtag'] + '/' + Info['f_string'] + '/')


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

import netCDF4 as nc

fn = f_dir + 'ocean_bry.nc'
#fn = '/Users/PM5/Documents/roms/forcing/ptx_2006/ocean_bry_1.nc'
ds = nc.Dataset(fn)
#zfun.ncd(ds)
f_ss = ds.variables['temp_south'][0,:,:].squeeze()
ds.close()

fn = r_dir + 'ocean_his_0025.nc'
ds = nc.Dataset(fn)
#zfun.ncd(ds)
r_ss = ds.variables['temp'][0,:,0,:].squeeze()
ds.close()

import matplotlib.pyplot as plt

plt.close()

fig = plt.figure(figsize=(16, 10))
cmap = plt.get_cmap(name='rainbow')

ax = fig.add_subplot(121)
cs = ax.pcolormesh(f_ss, cmap=cmap, vmin=0, vmax=10)
fig.colorbar(cs)

ax = fig.add_subplot(122)
cs = ax.pcolormesh(r_ss, cmap=cmap, vmin=0, vmax=10)
fig.colorbar(cs)

plt.show()

f0 = f_ss[0,:]
r0 = r_ss[0,:]
print(f0[75:85])
print(r0[75:85])



