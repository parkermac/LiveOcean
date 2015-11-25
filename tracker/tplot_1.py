"""
Plot results of tracker.

Designed for multiple releases and incorporating a mooring record.
"""

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun
import matfun
import matplotlib.pyplot as plt
import numpy as np

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'tracks_2014_CERF/'
 
import pickle
jdf, cr, G = pickle.load( open( indir + 'starters_2014.p', 'rb' ) )

NT, NP = jdf['lon'].shape

lonp = G['lon_psi']
latp = G['lat_psi']
aa_jdf = [-126, -123.5, 47, 49]
aa_cr = [-125, -122.5, 45, 47]
aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]   
depth_levs = [100, 200, 500, 1000, 2000, 3000]
    
fn_coast = Ldir['data'] + 'coast/pnw_coast_combined.mat'
cmat = matfun.loadmat(fn_coast)
    
# PLOTTING

cs_divider = -.4
csd_text = str(int(np.abs(100.*cs_divider)))
         
plt.close()
fig = plt.figure(figsize=(16,8))

# MAP OF TRACKS

ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

riv_list = ['jdf', 'cr']
for riv in riv_list:
    if riv == 'jdf':
        pp = jdf
        ax = ax1
        aa = aa_jdf
        title_text = 'Strait of Juan de Fuca'
    elif riv == 'cr':
        pp = cr
        ax = ax2
        aa = aa_cr   
        title_text = 'Columbia River'

    ax.contour(G['lon_rho'], G['lat_rho'], G['h'], depth_levs, colors='g')        
    ax.plot(cmat['lon'],cmat['lat'], '-k') # coastline       
    ax.axis(aa)
    zfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    if riv == 'jdf':
        ax.text(.85, .25, 'Depth Contours', horizontalalignment='center', transform=ax.transAxes, color='g')
        dd = 1
        for d in depth_levs:
            ax.text(.85, .25 - .03*dd, str(d), horizontalalignment='center', transform=ax.transAxes, color='g')
            dd += 1
    
    if False:  
        # plot starting points (assumes these are all reverse runs)
        # Color by end cs.
        ii = 0  
        for cs in pp['cs1'][0,:]:
            if cs > cs_divider:
                ax.plot(pp['lon'][:,ii],pp['lat'][:,ii],'or',markersize=5, alpha = .4, markeredgecolor='r')
                pass
            else:
                ax.plot(pp['lon'][:,ii],pp['lat'][:,ii],'ob',markersize=5, alpha = .4, markeredgecolor='b')
            ii += 1
        if riv == 'cr':
            ax.text(.95,.9, 'End depth above ' + csd_text + '%',
                horizontalalignment='right', transform=ax.transAxes,
                color='r', fontsize=16)
            ax.text(.95,.8, 'End depth below ' + csd_text + '%',
                horizontalalignment='right', transform=ax.transAxes,
                color='b', fontsize=16)
    else:  
        # plot starting points (assumes these are all reverse runs)
        # Color by starting cs
        csf = pp['cs'].flatten()
        lonf = pp['lon'].flatten()
        latf = pp['lat'].flatten()
        mask = csf <= -.7
        ax.plot(lonf[mask],latf[mask],'ob',markersize=5, alpha = .4, markeredgecolor='b')
        mask = csf >= -.3
        ax.plot(lonf[mask],latf[mask],'or',markersize=5, alpha = .4, markeredgecolor='r')
        mask = (csf < -.3) & (csf > -.7)
        ax.plot(lonf[mask],latf[mask],'og',markersize=5, alpha = .2, markeredgecolor='g')
        
        if riv == 'cr':
            ax.text(.95,.9, 'Start depth in top 30%',
                horizontalalignment='right', transform=ax.transAxes,
                color='r', fontsize=16)
            ax.text(.95,.85, 'Start in between',
                horizontalalignment='right', transform=ax.transAxes,
                color='g', fontsize=16)
            ax.text(.95,.8, 'Start depth in bottom 30%',
                horizontalalignment='right', transform=ax.transAxes,
                color='b', fontsize=16)
        
    # ending point
    ax.plot(pp['lon1'][0,0],pp['lat1'][0,0],'y*',markersize=20)
    ax.set_title(title_text)

plt.show()

