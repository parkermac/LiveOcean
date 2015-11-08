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
P, G, S, PLdir = pickle.load( open( indir + inname, 'rb' ) )

NT, NP = P['lon'].shape

lonp = G['lon_psi']
latp = G['lat_psi']
if 'jdf' in inname:
    aa = [-126, -123.5, 47.5, 49.]
elif 'cr' in inname:
    aa = [-125, -122.5, 45, 47]
else:
    aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]   
depth_levs = [100, 200, 500, 1000, 2000, 3000]
    
# get coastline
cmat = matfun.loadmat(fn_coast)
    
# PLOTTING
         
plt.close()
fig = plt.figure(figsize=(16,8))

# MAP OF TRACKS
   
ax = fig.add_subplot(121)
ax.contour(G['lon_rho'], G['lat_rho'], G['h'], depth_levs, colors='g')        
ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5) # coastline       
ax.axis(aa)
zfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# plot tracks
cs_divider = -.6

ii = 0
if PLdir['dir_tag'] == 'forward':
    time_index_for_color = 0
elif PLdir['dir_tag'] == 'reverse':
    time_index_for_color = NT-1    
for cs in P['cs'][time_index_for_color,:]:
    if cs < cs_divider:
        ax.plot(P['lon'][:, ii],P['lat'][:, ii],'-b')
    else:
        ax.plot(P['lon'][:, ii],P['lat'][:, ii],'-r')
    ii += 1
    
# starting points
ii = 0
for cs in P['cs'][time_index_for_color,:]:
    if cs < cs_divider:
        ax.plot(P['lon'][0,ii],P['lat'][0,ii],'ob',markersize=15)
    else:
        ax.plot(P['lon'][0,ii],P['lat'][0,ii],'or',markersize=8)
    ii += 1
    
# ending points
ax.plot(P['lon'][-1,:],P['lat'][-1,:],'y*',markersize=20)
  
ax.set_title(PLdir['dir_tag'].title())

# TIME SERIES
tdays = (P['ot'] - P['ot'][0])/86400.

ax = fig.add_subplot(3,2,2)
ii = 0
if PLdir['dir_tag'] == 'forward':
    time_index_for_color = 0
elif PLdir['dir_tag'] == 'reverse':
    time_index_for_color = NT-1    
for cs in P['cs'][time_index_for_color,:]:
    if cs < cs_divider:
        ax.plot(tdays, P['salt'][:,ii],'-b')
    else:
        ax.plot(tdays, P['salt'][:,ii],'-r')
    ii += 1
ax.set_ylabel('Salinity')

ax = fig.add_subplot(3,2,4)
ii = 0
if PLdir['dir_tag'] == 'forward':
    time_index_for_color = 0
elif PLdir['dir_tag'] == 'reverse':
    time_index_for_color = NT-1    
for cs in P['cs'][time_index_for_color,:]:
    if cs < cs_divider:
        ax.plot(tdays, P['temp'][:,ii],'-b')
    else:
        ax.plot(tdays, P['temp'][:,ii],'-r')
    ii += 1
ax.set_ylabel('Temperature $^{\circ}C$')

ax = fig.add_subplot(3,2,6)
ii = 0
if PLdir['dir_tag'] == 'forward':
    time_index_for_color = 0
elif PLdir['dir_tag'] == 'reverse':
    time_index_for_color = NT-1    
for cs in P['cs'][time_index_for_color,:]:
    if cs < cs_divider:
        ax.plot(tdays, P['z'][:,ii],'-b')
    else:
        ax.plot(tdays, P['z'][:,ii],'-r')
    ii += 1
ax.set_xlabel('Days')
ax.set_ylabel('Z (m)')
plt.show()