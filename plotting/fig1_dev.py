"""
Code to test extracting a three-year salinity time series from
Ecology casts, to use in fig1.py.
"""

import pandas as pd
import matplotlib.pyplot as plt


plt.close('all')
fig = plt.figure(figsize=(13,8))

ax = fig.add_subplot(211)

clist = list('bgr')

S = pd.DataFrame(columns=['S_bot', 'S_top'])
indir = '/Users/pm7/Documents/ptools_data/ecology/'
counter = 0
for year in [2017, 2018, 2019]:
    
    Casts = pd.read_pickle(indir + 'Casts_' + str(year) + '.p')
    casts = Casts[Casts['Station']=='HCB010']
    casts = casts.set_index('Date')
    casts = casts[['Salinity', 'Z']]
    # identify a single cast by its DATE
    calldates = casts.index
    dates = calldates.unique() # a short list of unique dates (1 per cast)
    for dd in dates:
        # NOTE the brackets around [dd] keep the result as a DataFrame even if
        # we are only pulling out a single row.
        ca = casts.loc[[dd],:]
        S.loc[dd,'S_bot'] = ca.iloc[0].Salinity
        S.loc[dd,'S_top'] = ca.iloc[-1].Salinity
        ca.plot(ax=ax, x='Salinity', y='Z', color=clist[counter], legend=False)
    counter += 1
ax = fig.add_subplot(212)
S.plot(ax=ax, style='*')
        
plt.show()
    
    