"""
Plot results of tracker.
"""

# setup
import os; import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path: sys.path.append(alp)
import Lfun; reload(Lfun)
import zfun; reload(zfun) # plotting functions
import matfun; reload(matfun) # functions for working with mat files
import matplotlib.pyplot as plt

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'tracks/'
fn_coast = Ldir['data'] + 'coast/pnw_coast_combined.mat'

# choose the type of plot to make
print '\n%s\n' % '** Choose mooring file to plot **'
m_list_raw = os.listdir(indir)
m_list = []
for m in m_list_raw:
    if m[-2:] == '.p':
        m_list.append(m)
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print str(npt) + ': ' + m_list[npt]
my_npt = int(raw_input('-- Input number -- '))
inname = m_dict[my_npt]
  
import cPickle as pickle
Plon, Plat, Pcs, G, PLdir = pickle.load( open( indir + inname, 'rb' ) )

lonp = G['lon_psi']
latp = G['lat_psi']
aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]   
depth_levs = [100, 200, 500, 1000, 2000, 3000]
    
# get coastline
cmat = matfun.loadmat(fn_coast)
    
# PLOTTING       
#plt.close()
fig = plt.figure(figsize=(8,16))

# 1. surface salinity    
ax = fig.add_subplot(111)
ax.contour(G['lon_rho'], G['lat_rho'], G['h'], depth_levs, colors='g')        
# add coastline
if len(fn_coast) != 0:
    ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)        
# extras
ax.axis(aa)
zfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

if 'reverse' in inname:
    ax.plot(Plon[:, Pcs[-1,:] <= - 0.5],Plat[:, Pcs[-1,:] < - 0.4],'-b')
    ax.plot(Plon[:, Pcs[-1,:] >= - 0.3],Plat[:, Pcs[-1,:] >= - 0.4],'-r')
else:
    ax.plot(Plon[:, Pcs[0,:] <= - 0.5],Plat[:, Pcs[0,:] < - 0.4],'-b')
    ax.plot(Plon[:, Pcs[0,:] >= - 0.3],Plat[:, Pcs[0,:] >= - 0.4],'-r')
ax.plot(Plon[0,:],Plat[0,:],'ok')

ax.set_title(PLdir['dir_tag'].title())

plt.show()