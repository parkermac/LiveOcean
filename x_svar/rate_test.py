"""
Code to test the calculation of ds/dt using the disgnostics, averages and history files.
"""

# imports
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta

import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
import zrfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
ex_name = 'lo8da'
Ldir['gtagex'] = Ldir['gtag'] + '_' + ex_name 

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

# names of input files
dtf = datetime(2018,6,10)
f_string = 'f' + dtf.strftime('%Y.%m.%d')
dir0 = Ldir['roms'] + '/output/' + Ldir['gtagex'] + '/' + f_string + '/'

n = 2 # avg, dia file number to look at
N0 = ('0000' + str(n))[-4:]
N1 = ('0000' + str(n+1))[-4:]
h0 = dir0 + 'ocean_his_'+N0+'.nc'
h1 = dir0 + 'ocean_his_'+N1+'.nc'
a = dir0 + 'ocean_avg_'+N0+'.nc'
d = dir0 + 'ocean_dia_'+N0+'.nc'

# initialize Datasets
dsh0 = nc.Dataset(h0)
dsh1 = nc.Dataset(h1)
dsa = nc.Dataset(a)
dsd = nc.Dataset(d)

# get fields (one vertical column)
s0 = dsh0['salt'][0,:,10,10]
s1 = dsh1['salt'][0,:,10,10]
r = dsd['salt_rate'][0,:,10,10]
s = dsa['salt'][0,:,10,10]
e0 = dsh0['zeta'][0,10,10]
e1 = dsh1['zeta'][0,10,10]
e = dsa['zeta'][0,10,10]
h = dsh0['h'][10,10]

S = zrfun.get_basic_info(h1, only_S=True)
zr0, zw0 = zrfun.get_z(h, e0, S)
zr1, zw1 = zrfun.get_z(h, e1, S)
zr, zw = zrfun.get_z(h, e, S)
dz0 = np.diff(zw0)
dz1 = np.diff(zw1)
dz = np.diff(zw)

# compute dsdt

dsdt = (s1 - s0)/3600
dsdt = zfun.fillit(dsdt)

dsdt_alt = r - s*(dz1 - dz0)/(3600*dz)
dsdt_alt = zfun.fillit(dsdt_alt)

err = dsdt_alt - dsdt

print('Error = %g out of %g' % (err.std(), dsdt.std()))

