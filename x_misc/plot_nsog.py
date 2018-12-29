"""
Plot aimed at understanding the Northern SoG boundary condition.

Pulls in data from a variety of sources.
"""

# setup
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
import pickle
import pandas as pd

import os
import sys
pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zfun

# get the nsog TEF section
pth = os.path.abspath(Ldir['LO'] + 'x_tef/')
if pth not in sys.path:
    sys.path.append(pth)
import tef_fun
indir = Ldir['LOo'] + 'tef/cas4_v2_lo6biom_2017.01.01_2017.12.31/'
fn = indir + 'sog4.p'
Qi, Si, Fi, qnet_lp, fnet_lp, td = tef_fun.tef_integrals(fn)
tef_dt = []
for tt in td:
    tef_dt.append(datetime(2017,1,1) + timedelta(days=tt - 0.5))
tef_ser = pd.Series(index=tef_dt, data=qnet_lp)

# get the Campbell River tide extraction
pth = os.path.abspath(Ldir['parent'] + 'ptools/tide_obs_mod/')
if pth not in sys.path:
    sys.path.append(pth)
import obsfun
indir = Ldir['parent'] + 'ptools_output/tide/'
noaa_sn_dict, dfo_sn_dict, sn_dict = obsfun.get_sn_dicts()
# load data
year  = 2017
name = 'Campbell River'
sn = sn_dict[name]
mod_dir = indir + 'mod_data/cas4_v2_lo6biom/'
fn = mod_dir + 'tide_' + str(sn) + '_' + str(year) + '.p'
tide = pd.read_pickle(fn)
eta = np.array(tide['eta'].tolist())
eta_lp = zfun.filt_godin(eta)
tide['eta_lp'] = eta_lp
tlp = tide['eta_lp']
tlpd = tlp.resample('1D').mean()
tlpd = tlpd.tz_localize(None)

# get the hycom field from nsog
indir = Ldir['LOo'] + 'misc/'
fn = indir + 'zeta_df.p'
zdf = pd.read_pickle(fn)

# make a DataFrame that has all fields on the same time axis

df = pd.DataFrame(index=pd.date_range(start='1/1/2017', end='1/1/2018'))
df['ROMS Campbell River LP SSH'] = tlpd
df['HYCOM N SoG SSH'] = zdf['z_sog']
df['HYCOM JdF SSH'] = zdf['z_jdf']
df['SSH HYCOM-ROMS'] = df['HYCOM N SoG SSH'] - df['ROMS Campbell River LP SSH']
df['TEF N SoG Q'] = tef_ser

# plotting
plt.close('all')

fig = plt.figure(figsize=(14,8))

ax = fig.add_subplot(221)
df.plot(y=['ROMS Campbell River LP SSH', 'HYCOM N SoG SSH', 'HYCOM JdF SSH'], ax=ax)
ax.set_xlim(datetime(2017,1,1), datetime(2018,1,1))

ax = fig.add_subplot(223)
df.plot(y='TEF N SoG Q', ax=ax)
ax.set_xlim(datetime(2017,1,1), datetime(2018,1,1))

ax = fig.add_subplot(222)
df.plot(x='SSH HYCOM-ROMS', y='TEF N SoG Q', style='*k', ax=ax, grid=True)

ax = fig.add_subplot(224)
df.plot(x='HYCOM JdF SSH', y='HYCOM N SoG SSH', style='ob', ax=ax, grid=True)
ax.plot([0, .3], [0, .3], '-k')

plt.show()

