"""
Compares two history files.
"""

# setup
import os; import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(alp)
import zfun; reload(zfun)


f_string = 'f2013.01.02'
nhiss = '0025'
fn1 = (Ldir['roms'] + 'output/' + Ldir['gtag'] + '/'
    + f_string + '/ocean_his_' + nhiss + '.nc')
fn2 = (Ldir['roms'] + 'output/' + Ldir['gtag'] + '_v0/'
    + f_string + '/ocean_his_' + nhiss + '.nc')
    
# GET DATA
G, S, T = zfun.get_basic_info(fn1)
import netCDF4 as nc   
ds = nc.Dataset(fn1,'r')
h = G['h']
salt1 = ds.variables['salt'][0, -1, :, :].squeeze()
temp1 = ds.variables['temp'][0, -1, :, :].squeeze()
u1 = ds.variables['u'][0, -1, :, :].squeeze()
v1 = ds.variables['v'][0, -1, :, :].squeeze()  
ds.close()     
lonp = G['lon_psi']
latp = G['lat_psi']
aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]

ds = nc.Dataset(fn2,'r')
salt2 = ds.variables['salt'][0, -1, :, :].squeeze()
temp2 = ds.variables['temp'][0, -1, :, :].squeeze()
u2 = ds.variables['u'][0, -1, :, :].squeeze()
v2 = ds.variables['v'][0, -1, :, :].squeeze()  
ds.close()      
        
# PLOTTING       
import matplotlib.pyplot as plt
plt.close()
fig = plt.figure(figsize=(14, 8))

# 1. surface salinity    
ax = fig.add_subplot(121)
cmap = plt.get_cmap(name='jet')    
cs = ax.pcolormesh(lonp, latp,
    salt1[1:-1,1:-1] - salt2[1:-1,1:-1],
    cmap = cmap)
ax.axis(aa)
zfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Surface Salinity Difference')
fig.colorbar(cs)
# put info on plot
ax.text(.95, .1, T['tm'].strftime('%Y-%m-%d'),
    horizontalalignment='right', transform=ax.transAxes)
ax.text(.95, .05, T['tm'].strftime('%H:%M'),
    horizontalalignment='right', transform=ax.transAxes)

# 2. Temperature
ax = fig.add_subplot(122)
cmap = plt.get_cmap(name='jet')    
cs = ax.pcolormesh(lonp, latp,
    temp1[1:-1,1:-1] - temp2[1:-1,1:-1],
    cmap = cmap)
ax.axis(aa)
zfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Surface Temperature Difference')
fig.colorbar(cs)

plt.show()   
