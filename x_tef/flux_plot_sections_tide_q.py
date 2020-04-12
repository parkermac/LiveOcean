"""
Plot the mean of tidal energy flux and volume transport at all
TEF sections.

"""

# imports
import matplotlib.pyplot as plt
import pickle
import netCDF4 as nc
import pandas as pd
import numpy as np

import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
Ldir = Lfun.Lstart(gridname, tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + ex_name

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

import tef_fun
import flux_fun
from importlib import reload
reload(flux_fun)

# User choices
year = 2018
year_str = str(year)
testing = True

# input
run_name = Ldir['gtagex']+'_'+year_str+'.01.01_'+year_str+'.12.31'
indir00 = Ldir['LOo'] + 'tef/'
indir0 = indir00 + run_name + '/'
indir = indir0 + 'bulk/'

# output
out_fn = indir00 + 'misc_figs_' + Ldir['gridname'] + '/Tide_and_Qnet_' + year_str + '.png'

# colors
clist = flux_fun.clist

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()
sect_list = list(sect_df.index)

# if testing:
#     sect_list = [sect_list[0]]
    
df = pd.DataFrame(index=sect_list)
for sn in sect_list:
    bulk = pickle.load(open(indir + sn + '.p', 'rb'))
    sx0, sx1, sy0, sy1, landward = sect_df.loc[sn,:]
    sx = (sx0+sx1)/2; sy = (sy0+sy1)/2
    if (sx0==sx1) and (sy0!=sy1):
        sdir = 'NS'
    elif (sx0!=sx1) and (sy0==sy1):
        sdir = 'EW'
    df.loc[sn,'lon'] = sx
    df.loc[sn,'lat'] = sy
    df.loc[sn,'sdir'] = sdir
    df.loc[sn,'landward'] = landward
    df.loc[sn,'F'] = bulk['fnet_lp'].mean()/1e6 # MW
    df.loc[sn,'Q'] = bulk['qnet_lp'].mean() # 1000 m3/s
    
# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(10,13))
fs = 18
cf = 'purple'
cq = 'g'

# axis limits
x0 = -125.5; x1 = -122; y0 = 47; y1 = 50.5 # Salish Sea
x00 = -123.3; x11 = -122.2; y00 = 47; y11 = 48.5 # Puget Sound
aaS = [x0, x1, y0, y1]
aaP = [x00, x11, y00, y11]

ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

for sn in df.index:
    sx = df.loc[sn,'lon']
    sy = df.loc[sn,'lat']
    sdir = df.loc[sn,'sdir']
    landward = df.loc[sn,'landward']
    # landward is the sign to multiply by to know what the direction
    # of a positive flux is.  So if sdir = 'NS' and landward = 1,
    # then positive flux is to the East
    F = df.loc[sn,'F']
    Q = df.loc[sn,'Q']
    logF = np.log10(np.abs(F) + 1)
    logQ = np.log10(np.abs(Q) + 1)
    sgnF = np.sign(landward*F)
    sgnQ = np.sign(landward*Q)
    scl = 20 # a vector of length 1 will be 1/scl of the length of the y-axis
    if sdir == 'EW':
        vf0 = 0; vf1 = sgnF
        vq0 = 0; vq1 = sgnQ
    elif sdir == 'NS':
        vf0 = sgnF; vf1 = 0
        vq0 = sgnQ; vq1 = 0
        
    do_text = True
    ax1.quiver(sx,sy, vf0, vf1, scale=scl, scale_units='height',
        headwidth=0, headlength=0, linewidths=1, color=cf)
    ax1.plot(sx,sy, 'ok', markersize=10*logF, markerfacecolor='None', markeredgecolor=cf)
    if do_text:
        ax1.text(sx, sy+.04, str(int(np.abs(F))), ha='center', va='center', size=.5*fs,
        weight='bold', color=cf, alpha=.8)
        
    ax2.quiver(sx,sy, vf0, vf1, scale=scl, scale_units='height',
        headwidth=0, headlength=0, linewidths=1, color=cf)
    ax2.plot(sx,sy, 'ok', markersize=10*logF, markerfacecolor='None', markeredgecolor=cf)
    if do_text and sx > x00 and sy < y11:
        ax2.text(sx, sy+.02, str(int(np.abs(F))), ha='center', va='center', size=.5*fs,
        weight='bold', color=cf, alpha=.8)

    ax3.quiver(sx,sy, vq0, vq1, scale=scl, scale_units='height',
        headwidth=0, headlength=0, linewidths=1, color=cq)
    ax3.plot(sx,sy, 'ok', markersize=10*logQ, markerfacecolor='None', markeredgecolor=cq)
    if do_text:
        ax3.text(sx, sy+.04, str(int(np.abs(Q))), ha='center', va='center', size=.5*fs,
        weight='bold', color=cq, alpha=.8)

    ax4.quiver(sx,sy, vq0, vq1, scale=scl, scale_units='height',
        headwidth=0, headlength=0, linewidths=1, color=cq)
    ax4.plot(sx,sy, 'ok', markersize=10*logQ, markerfacecolor='None', markeredgecolor=cq)
    if do_text and sx > x00 and sy < y11:
        ax4.text(sx, sy+.02, str(int(np.abs(Q))), ha='center', va='center', size=.5*fs,
        weight='bold', color=cq, alpha=.8)

pfun.add_coast(ax1, color='gray')
ax1.axis(aaS)
pfun.dar(ax1)
ax1.set_xticklabels([])
ax1.tick_params(labelsize=.8*fs)
ax1.text(.95,.9,'(a)', size=fs, transform=ax1.transAxes, ha='right', weight='bold')
ax1.set_ylabel('Latitude', size=.8*fs)
ax1.set_xticks([-125, -124, -123, -122])
for FF in [1, 100, 10000]:
    ax1.plot(-125, 47.5, 'o', markersize=10*np.log10(np.abs(FF) + 1),
    markerfacecolor='None', markeredgecolor=cf)
ax1.text(-125.4, 47.15, 'Tidal Energy Flux\n(1, 100, 10000) [MW]',
    ha='left', va='center', size=.6*fs)

pfun.add_coast(ax2, color='gray')
ax2.axis(aaP)
pfun.dar(ax2)
ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.tick_params(labelsize=.8*fs)
ax2.text(.95,.9,'(b)', size=fs, transform=ax2.transAxes, ha='right', weight='bold')
ax2.set_xticks([-123, -122.5])

pfun.add_coast(ax3, color='gray')
ax3.axis(aaS)
pfun.dar(ax3)
ax3.tick_params(labelsize=.8*fs)
ax3.text(.95,.9,'(c)', size=fs, transform=ax3.transAxes, ha='right', weight='bold')
ax3.set_xlabel('Longitude', size=.8*fs)
ax3.set_ylabel('Latitude', size=.8*fs)
for QQ in [1, 100, 10000]:
    ax3.plot(-125, 47.5, 'o', markersize=10*np.log10(np.abs(QQ) + 1),
    markerfacecolor='None', markeredgecolor=cq)
ax3.text(-125.4, 47.15, 'Volume Transport\n(1, 100, 10000) [$m^{3}s^{-1}$]',
    ha='left', va='center', size=.6*fs)
ax3.set_xticks([-125, -124, -123, -122])
ax3.set_xticklabels([-125, -124, -123, -122])

pfun.add_coast(ax4, color='gray')
ax4.axis(aaP)
pfun.dar(ax4)
ax4.set_yticklabels([])
ax4.tick_params(labelsize=.8*fs)
ax4.text(.95,.9,'(d)', size=fs, transform=ax4.transAxes, ha='right', weight='bold')
ax4.set_xlabel('Longitude', size=.8*fs)
ax4.set_xticks([-123, -122.5])
ax4.set_xticklabels([-123, -122.5])


        
fig.tight_layout()
plt.savefig(out_fn)
plt.show()


    

    