"""
Code to see what the observed 3-year trend in salinity is.

"""

# imports

#%% setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart()
import zfun, zrfun
import pfun

from datetime import datetime, timedelta
import numpy as np
import netCDF4 as nc
import pickle
import pandas as pd
import matplotlib.pyplot as plt

dir0 = Ldir['parent'] + 'ptools_data/ecology/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(14,7))

dt0 = datetime(2017,1,1,0)
dt1 = datetime(2020,1,1,0)

ax = fig.add_subplot(111)

# add observations
sta_list = [item for item in sta_df.index if 'SJF' in item]
for station in sta_list:
    S = pd.DataFrame(columns=['S_bot', 'S_top'])
    indir = '/Users/pm7/Documents/ptools_data/ecology/'
    for year in [2017, 2018, 2019]:
        Casts = pd.read_pickle(indir + 'Casts_' + str(year) + '.p')
        casts = Casts[Casts['Station']==station]
        casts = casts.set_index('Date')
        casts = casts[['Salinity', 'Z']]
        # identify a single cast by its DATE
        calldates = casts.index
        dates = calldates.unique() # a short list of unique dates (1 per cast)
        for dd in dates:
            # NOTE the brackets around [dd] keep the result as a DataFrame even if
            # we are only pulling out a single row.
            ca = casts.loc[[dd],:]
            ss = ca[ca['Z']==-50].Salinity.to_numpy()
            if len(ss) == 1:
                S.loc[dd,'S_bot'] = ss[0]
            else:
                S.loc[dd,'S_bot'] = np.nan
            S.loc[dd,'S_top'] = ca.iloc[-1].Salinity
    if len(S) > 0:
        S['S_bot'].plot(ax=ax, style='-b*', grid=True, legend=False)
        S['S_top'].plot(ax=ax, style='-r*', grid=True, legend=False)

plt.show()



