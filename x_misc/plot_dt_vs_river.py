"""
Code to plot the DT(day) of LiveOcean, comparing it
to selected river flow.
"""

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun

import pandas as pd
import matplotlib.pyplot as plt

Ldir = Lfun.Lstart()

fn_dt = Ldir['LOo'] + 'misc/dt_log_cas4_v2_lo6biom.txt'
df_dt = pd.read_csv(fn_dt,parse_dates=True,index_col='Date')

fn_r = Ldir['LOo'] + 'river/cas4_v2_2017.01.01_2018.12.31.p'
df_r = pd.read_pickle(fn_r)

# also get superplot info to look at tides
ffn = Ldir['LOo'] + 'superplot/forcing_cas4_v2_lo6biom_2017.p'
df_super = pd.read_pickle(ffn)

plt.close('all')


df_comb = pd.DataFrame(index=df_r.index, columns=['DT', 'fraser', 'columbia'])
df_comb['fraser'] = df_r['fraser']
df_comb['columbia'] = df_r['columbia']
df_comb['DT'] = df_dt['DT']

# merge in super tides
df_comb = df_comb.loc[df_super.index,:]

df_comb['RMS Tide Height (m)'] = df_super['RMS Tide Height (m)']

fig = plt.figure(figsize=(14,8))

t0 = df_comb.index[0]
t1 = df_comb.index[-1]

ax = fig.add_subplot(221)
df_comb['DT'].plot(ax=ax, xlim=(t0,t1))
ax.set_ylabel('DT (sec)')
ax.set_title('cas4_v2_lo6biom')

ax = fig.add_subplot(223)
df_comb['fraser'].plot(ax=ax, xlim=(t0,t1), style='-r', label='fraser')
df_comb['columbia'].plot(ax=ax, xlim=(t0,t1), style='-b', label='columbia')
ax.legend()

ax = fig.add_subplot(222)
df_comb.plot(x='RMS Tide Height (m)', y='DT', style='om', markersize=10,
    alpha=.3, ax=ax, legend=False)
ax.set_ylabel('DT (sec)')

ax = fig.add_subplot(224)
df_comb.plot(x='fraser', y='DT', style='or', markersize=10,
    alpha=.3, ax=ax, legend=False)
df_comb.plot(x='columbia', y='DT', style='ob', markersize=5,
    alpha=.3, ax=ax, legend=False)
ax.set_ylabel('DT')
ax.set_xlabel('River Flow (m3/s)')

plt.show()