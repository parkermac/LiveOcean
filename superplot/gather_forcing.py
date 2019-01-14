"""
Code to gather time series related to LiveOcean forcing.
"""

import os
import sys
pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zrfun
import zfun

import numpy as np
import pandas as pd
import pickle
import netCDF4 as nc
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

# define where to put the output
outdir = Ldir['LOo'] + 'superplot/'
Lfun.make_dir(outdir)

gtag = 'cas4_v2'
gtagex = gtag + '_lo6biom'
year = '2017'

# River forcing
# created by x_river/extract_rivers.py
fnr = gtag + '_' + year + '.01.01_' + year + '.12.31.p'
fn = Ldir['LOo'] + 'river/' + fnr
riv_df0 = pd.read_pickle(fn)
# keep only a subset of the rivers
riv_df = riv_df0.loc[:, ['columbia', 'fraser', 'skagit']]

# Tide
pth = os.path.abspath(Ldir['parent'] + 'ptools/tide_obs_mod/')
if pth not in sys.path:
    sys.path.append(pth)
import obsfun
noaa_sn_dict, dfo_sn_dict, sn_dict = obsfun.get_sn_dicts()
t_indir = Ldir['parent'] + 'ptools_output/tide/mod_data/' + gtagex + '/'
t_fn = t_indir + 'tide_' + str(sn_dict['Seattle']) + '_' + year + '.p'
tide_df = pickle.load(open(t_fn, 'rb'))
# remove the timezone
tide_df = tide_df.tz_localize(None)
# make a dataframe of spring-neap conditions
eta = tide_df['eta'].values
eta_rms = np.sqrt(zfun.filt_godin(eta**2))
tide_df['eta_rms'] = eta_rms
# subsample to daily
tide_daily_df = tide_df.loc[::24, 'eta_rms']

# Wind
fnw = Ldir['LOo'] + 'moor/' + gtagex + '_2017.01.01_2018.11.29/NANOOS_ChaBa_Buoy_hourly.nc'
moor_ds = nc.Dataset(fnw)
ot = moor_ds['ocean_time'][:]
svstr = moor_ds['svstr'][:]
svstr_lp = zfun.filt_AB8d(svstr)
wind_time = []
for tt in ot:
    wind_time.append(Lfun.modtime_to_datetime(tt))
wind_df = pd.DataFrame(index=wind_time, columns=['svstr', 'svstr_lp'])
wind_df['svstr'] = svstr
wind_df['svstr_lp'] = svstr_lp
wind_df = wind_df.resample('D').mean()

# combine daily fields into a single DataFrame
comb_df = pd.DataFrame(index=pd.date_range(start='1/1/2017', end='12/31/2017'))
comb_df['RMS Tide Height (m)'] = tide_daily_df
comb_df['8-day NS Wind Stress (Pa)'] = wind_df['svstr_lp']
comb_df['Columbia R. Flow (1000 m3/s)'] = riv_df['columbia']/1000
comb_df['Fraser R. Flow (1000 m3/s)'] = riv_df['fraser']/1000
comb_df['Skagit R. Flow (1000 m3/s)'] = riv_df['skagit']/1000

# save for later use
out_fn = outdir + 'forcing_' + gtagex + '_' + year + '.p'
print('Saving ' + out_fn)
comb_df.to_pickle(out_fn)

# plotting
plt.close('all')
comb_df.plot(subplots=True, xlim=(comb_df.index[0], comb_df.index[-1]), grid=True)
plt.show()
