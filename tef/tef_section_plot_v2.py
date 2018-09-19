"""
Plot the mean of many TEF extractions on a thalweg section.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import pickle

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun

import tef_fun
from importlib import reload
reload(tef_fun)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

Ldir = Lfun.Lstart()

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

indir = Ldir['LOo'] + 'tef/'
if False:
    print('\nSelect an Extraction to plot:\n')
    List_raw = os.listdir(indir)
    List_raw.sort()
    List = [item for item in List_raw]
    NL = len(List)
    Ldict = dict(zip(range(NL), List))
    for ii in range(NL):
        print(str(ii) + ': ' + List[ii])
    if True:
        my_ii = int(input('-- Input number: '))
    else:
        my_ii = 0 # for testing
    Litem = Ldict[my_ii]
else:
    Litem = 'cas4_v1_lo6biom_2017.01.01_2017.12.31'
print('\nProcessing ' + Litem + '\n')
LList_raw = os.listdir(indir + Litem)
LList_raw.sort()
LList = [item for item in LList_raw if ('.p' in item)]
Indir = indir + Litem + '/'

plt.close('all')

# override
LList = ['jdf1.p','jdf2.p','jdf3.p','jdf4.p',
    'ai1.p', 'ai2.p', 'ai3.p','ai4.p',
    'mb1.p','mb2.p','mb3.p','mb4.p','mb5.p',
    'tn1.p','tn2.p','tn3.p',
    'ss1.p','ss2.p','ss3.p']
NS = len(LList)
dd = np.nan * np.ones(NS)

qin = np.nan * np.ones(NS)
qout = np.nan * np.ones(NS)
qsin = np.nan * np.ones(NS)
qsout = np.nan * np.ones(NS)

qin_abs = np.nan * np.ones(NS)
qout_abs = np.nan * np.ones(NS)
qsin_abs = np.nan * np.ones(NS)
qsout_abs = np.nan * np.ones(NS)

sin = np.nan * np.ones(NS)
sout = np.nan * np.ones(NS)

xs = np.nan * np.ones(NS)
ys = np.nan * np.ones(NS)

counter = 0
for LL in LList:
    
    sect_name = LL.replace('.p','')
    print('\n** ' + sect_name + ' **')
    
    x0, x1, y0, y1, landward = sect_df.loc[sect_name,:]
    
    lon = (x0+x1)/2
    lat = (y0+y1)/2
    
    if counter == 0:
        lon0 = lon
        lat0 = lat
        
    xs[counter], ys[counter] = zfun.ll2xy(lon, lat, lon0, lat0)

    fn = Indir + LL
    
    Qi, Si, Fi, qnet_lp, fnet_lp, td = tef_fun.tef_integrals_v2(fn)
        
    ndg = 363 # number of good days
    
    
    qin[counter] = np.nansum(Qi[:,0]/1e3)/ndg
    qout[counter] = np.nansum(Qi[:,1]/1e3)/ndg
        
    qsin[counter] = np.nansum(Fi[:,0]/1e3)/ndg
    qsout[counter] = np.nansum(Fi[:,1]/1e3)/ndg
    
    qin_abs[counter] = np.nansum(np.abs(Qi[:,0]))/ndg
    qout_abs[counter] = np.nansum(np.abs(Qi[:,1]))/ndg
        
    qsin_abs[counter] = np.nansum(np.abs(Fi[:,0]))/ndg
    qsout_abs[counter] = np.nansum(np.abs(Fi[:,1]))/ndg
    
    sin[counter] = qsin_abs[counter]/qin_abs[counter]
    sout[counter] = qsout_abs[counter]/qout_abs[counter]
        
    counter += 1

# create a distance vector
dx = np.diff(xs)
dy = np.diff(ys)
dd = np.sqrt(dx**2 + dy**2)
dist = np.zeros(NS)
dist[1:] = np.cumsum(dd/1000)

# plotting

lw=2

fig = plt.figure(figsize=(13,8))

ax = fig.add_subplot(311)
ax.plot(dist,qin,'-or',linewidth=lw, label='Qin')
ax.plot(dist,-qout,'-ob',linewidth=lw, label='Qout')
ax.set_xlim(0,300)
ax.grid(True)
ax.set_ylabel('Transport (1000 m3/s)')
ax.legend()
counter = 0
for sn in LList:
    sn = sn.replace('.p','')
    sn = sn.upper()
    ax.text(dist[counter], qin[counter]+20, sn, rotation=45)
    counter += 1

ax = fig.add_subplot(312)
ax.plot(dist,qsin,'-*r',linewidth=lw, label='QSin')
ax.plot(dist,-qsout,'-*b',linewidth=lw, label='QSout')
ax.set_xlim(0,300)
ax.grid(True)
ax.set_ylabel('Salt Transport (psu 1000 m3/s)')
ax.legend()

ax = fig.add_subplot(313)
ax.plot(dist,sin,'-+r',linewidth=lw, label='Sin')
ax.plot(dist,sout,'-+b',linewidth=lw, label='Sout')
ax.set_xlim(0,300)
ax.grid(True)
ax.set_xlabel('Distance from Mouth (km)')
ax.set_ylabel('Salinity')
ax.legend()

plt.show()
    