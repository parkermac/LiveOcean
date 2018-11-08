"""
Code to look at the UBC extraction.
"""

import matplotlib.pyplot as plt
import netCDF4 as nc

fn = '/Users/pm7/Documents/LiveOcean_roms/output/cas4_v2_lo6biom/f2018.09.29/low_passed_UBC.nc'

ds = nc.Dataset(fn)

# for vn in ds.variables:
#     print(vn)

# plotting
plt.close('all')
fig = plt.figure(figsize=(12,8))

ax = fig.add_subplot(121)
a = ds['salt'][0,-1,:,:].squeeze()
x = ds['lon_rho'][:]
y = ds['lat_rho'][:]
cs = ax.pcolormesh(x,y,a)
fig.colorbar(cs)
ax.set_title('salt')

ax = fig.add_subplot(122)
a = ds['u'][0,-1,:,:].squeeze()
x = ds['lon_u'][:]
y = ds['lat_u'][:]
cs = ax.pcolormesh(x,y,a)
fig.colorbar(cs)
ax.set_title('u')

plt.show()