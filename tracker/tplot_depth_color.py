"""
Plot results of a particle tracking experiment.
"""

# setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
sys.path.append(os.path.abspath('../plotting'))
import pfun

import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np

Ldir = Lfun.Lstart()

# Choose an experiment to plot from.
indir0 = Ldir['LOo'] + 'tracks2/'
indir_list_raw = os.listdir(indir0)
indir_list = []
for d in indir_list_raw:
    if os.path.isdir(indir0 + d):
        indir_list.append(d)
indir_list.sort()
Npt = len(indir_list)#
print('\n%s\n' % '** Choose Experiment to plot **')
for npt in range(Npt):
    print(str(npt) + ': ' + indir_list[npt])
my_npt = input('-- Experiment number (return = 0) --')
if len(my_npt)==0:
    my_npt = 0
indir = indir_list[int(my_npt)] + '/'

# Choose a release from this experiment.
rel_list = [rel for rel in os.listdir(indir0 + indir) if 'release' in rel]
rel_list.sort()
Nrl = len(rel_list)
print('\n%s\n' % '** Choose Release file to plot **')
for nrl in range(Nrl):
    print(str(nrl) + ': ' + rel_list[nrl])
my_nrl = input('-- Release number (return = 0) -- ')
if len(my_nrl)==0:
    my_nrl = 0
rel = rel_list[int(my_nrl)]

# get Datasets
dsr = nc4.Dataset(indir0 + indir + rel)
dsg = nc4.Dataset(indir0 + indir + 'grid.nc')

NT, NP = dsr['lon'].shape

# get a list of datetimes
ot_vec = dsr['ot'][:]
dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]

# gather some fields, for convenience
lonp = dsg['lon_psi'][:]
latp = dsg['lat_psi'][:]
hh = dsg['h'][:]
maskr = dsg['mask_rho'][:] # 1=water, 0=land
#
u = dsr['u'][:]
v = dsr['v'][:]
w = dsr['w'][:]
salt = dsr['salt'][:]
temp = dsr['temp'][:]
lon = dsr['lon'][:]
lat = dsr['lat'][:]
z = dsr['z'][:]
cs = dsr['cs'][:]
zeta = dsr['zeta'][:]
h = dsr['h'][:]

# subsample output for plotting
npmax = 300 # max number of points to plot
step = np.max((int(NP/npmax),1))
u = u[:,::step]
v = v[:,::step]
w = w[:,::step]
salt = salt[:,::step]
temp = temp[:,::step]
lon = lon[:,::step]
lat = lat[:,::step]
z = z[:,::step]
cs = cs[:,::step]
zeta = zeta[:,::step]
h = h[:,::step]
NTS, NPS = lon.shape

# PLOTTING
#plt.close('all')
fig = plt.figure(figsize=(12,8))

# MAP
# set domain limits
if False:
    # plot full domain
    aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
else:
    # automatically plot region of particles, with padding
    pad = .02
    aa = [lon.min() - pad, lon.max() + pad,
    lat.min() - pad, lat.max() + pad]
ax = fig.add_subplot(121)
zm = -np.ma.masked_where(maskr==0, hh)
plt.pcolormesh(lonp, latp, zm[1:-1, 1:-1], vmin=-100, vmax=0,
    cmap='terrain', alpha=.25)
pfun.add_coast(ax)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title(indir.strip('/'))

# add the tracks (packed [time, particle])

# color end position by depth
z1 = z[-1,:]
x1 = lon[-1,:]
y1 = lat[-1,:]
sc_obj = ax.scatter(x1,y1, s=10, c=z1, cmap='rainbow')#,vmin=-20, vmax=5)
fig.colorbar(sc_obj)

# time series
td = (ot_vec - ot_vec[0])/86400
tv_list = ['z', 'salt', 'temp']
#tv_list = ['u', 'v', 'lon', 'lat']
ntv = len(tv_list)
for ii in range(ntv):
    tv = tv_list[ii]
    NC = 2
    ax = fig.add_subplot(ntv,NC, (ii+1)*NC)
    ax.plot(td, dsr[tv][:,::step])
    ax.text(.05, .05, tv, fontweight='bold', transform=ax.transAxes)
    if ii == ntv-1:
        ax.set_xlabel('Time (days)')

plt.show()

dsr.close()
dsg.close()

