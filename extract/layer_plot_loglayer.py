"""
Plots layer records with calculation of z at a specific height above
the bottom (assumes we are in the log layer).

Notes from the ROMS Forum from John Warner 7/13/2018:

"the routine set_vbc.F computes the bottom stresses
(in case you want to look at the code).

For QDRAG, it essentially uses (for u-direction)
bustr = rdrag2 * u(z=1) * |V|
For this option, as u(z=1) changes elevation, the roughness (rdrag2) stays the same.

Another option is to use UV_LOGDRAG
bustr = kappa^2 / (Ln(z/z0))^2 * u(z=1) * |V|
So when u(z=1) changes elevation, the model is assuming a log layer distribution
from the middle of the bottom cell to the sea floor,
this way the roughness scales with the elevation.

With your results, you could assume a log layer from the middle height
of the bottom cell to the sea floor. ROMS does not do anything fancy in
very shallow water where it may actually be resolving the vertical log layer profile,
unless you defined SPLINES. That can influence the shape of the velocity profile.
We typically do not use SPLINES in shallow water."

"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

# Rose - you will want to skip these lines and make your own indir
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'extract/'

# choose the mooring extraction to plot
print('\n%s\n' % '** Choose layer file to plot **')
m_list_raw = os.listdir(indir)
m_list_raw.sort()
m_list = [m for m in m_list_raw if (('.nc' in m) and ('layer_' in m))]
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
if True:
    my_npt = int(input('-- Input number -- '))
else:
    my_npt = 0 # for testing
moor_file = m_dict[my_npt]
fn = indir + moor_file

ds = nc.Dataset(fn)

# PLOTTING
plt.close('all')

fig = plt.figure(figsize=(16,8))

xp = ds['lon_psi'][:]
yp = ds['lat_psi'][:]
ot = ds['ocean_time'][:]

# log layer things
rho = 1026 # kg m-3

# this is the roughness length.  0.0008 is what I found to give
# a good fit meaning that when I used my bottom stresses and a log
# layer shape it gave values close to the original ones (see the scatter plot
# comparing uz_orig to uz_oldz)
z0 = 0.0007 # m (typical of mud/sand)
# for reference the model used a quadratic drag law with CD = 3e-3

# What I don't know is whether or not the model assumes a log layer profile
# when the layers are thin enough to resolve the log layer - e.g. in
# shallow water.

u = ds['u'][:]
v = ds['v'][:]
bustr = ds['bustr'][:]
bvstr = ds['bvstr'][:]
bstr = np.ma.sqrt(bustr**2 + bvstr**2)
u_star = np.ma.sqrt(bstr/rho)

zlay = ds['zlay'][:] # z position of velocities
zbot = -ds['h'][:] # z position of the bottom

zorig = zlay - zbot # m above the bottom of the velocities
zorig = np.ma.masked_where(u_star[0,:,:].mask, zorig)

# set the desired height above bottom at which you want to estimate velocities
znew = 0.1 * np.ones(zorig.shape)

def get_u_at_z(z, z0, u_star):
    # the classical log layer function
    uz = (u_star/0.41)* np.log(z/z0)
    return uz

# calculate three speeds:

# the original
uz_orig = np.ma.sqrt(u**2 + v**2)

# what our log layer predicts the speed should be at the
# original (varying) height above the bottom
uz_oldz = get_u_at_z(zorig, z0, u_star)

# what our log layer predicts the speed should be at the
# new (constant) height above the bottom
uz_newz = get_u_at_z(znew, z0, u_star)

# set up for plotting
F = dict() # field
S = dict() # series

F['uz_orig'] = uz_orig[0,:,:].squeeze()
F['uz_oldz'] = uz_oldz[0,:,:].squeeze()
F['uz_newz'] = uz_newz[0,:,:].squeeze()

ij = 50
S['uz_orig'] = uz_orig[:,ij,ij].squeeze()
S['uz_oldz'] = uz_oldz[:,ij,ij].squeeze()
S['uz_newz'] = uz_newz[:,ij,ij].squeeze()

v_list = ['uz_orig', 'uz_oldz', 'uz_newz']

days = (ot - ot[0])/86400.

NC = len(v_list)
count = 1
for vn in v_list:
    ax = fig.add_subplot(2,NC,count)
    cs = ax.pcolormesh(xp, yp, F[vn][1:-1, 1:-1],vmin=0, vmax=.4, cmap='rainbow')
    fig.colorbar(cs, ax=ax)
    ax.set_title(vn)
    count += 1

ax = fig.add_subplot(2,1,2)
for vn in v_list:
    ax.plot(days, S[vn],label=vn)
    ax.set_title(vn)
ax.legend()
plt.suptitle(fn)

fig2 = plt.figure(figsize=(8,8))
ax2 = fig2.add_subplot(111)
ax2.plot(uz_orig[0,:,:].flatten(),uz_oldz[0,:,:].flatten(),'.g',alpha=.2)
ax2.plot([0,1],[0,1],'-k')
ax2.grid(True)
ax2.set_xlabel('uz_orig')
ax2.set_ylabel('uz_oldz')

fig3 = plt.figure(figsize=(8,8))
ax3 = fig3.add_subplot(111)
ax3.plot(uz_orig[0,:,:].flatten(),uz_newz[0,:,:].flatten(),'.b',alpha=.2)
ax3.plot([0,1],[0,1],'-k')
ax3.grid(True)
ax3.set_xlabel('uz_orig')
ax3.set_ylabel('uz_newz')


plt.show()


