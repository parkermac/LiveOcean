"""
Plot Fraser River processed data.
"""

import os
import sys
pth = os.path.abspath('../../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import matplotlib.pyplot as plt

indir = Ldir['data'] + 'rivers/data_processed/'
qdf = pd.read_pickle(indir + 'fraser_flow.p')
qdfh = pd.read_pickle(indir + 'fraser_flow_historical.p')
qdfh2 = pd.read_pickle(indir + 'fraser_flow_historical_2.p')
tdf = pd.read_pickle(indir + 'fraser_temperature.p')

plt.close()
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(15,8), squeeze=False)

ax = axes[0, 0]
qdf.plot(ax=ax, title='Fraser River', style='-k')
qdfh.plot(ax=ax, style='-b')
qdfh2.plot(ax=ax, style='-r')
ax.set_ylabel('Flow (m3/s)')
ax.xaxis.set_ticklabels([])
ax.set_xlim(qdfh.index[0], qdf.index[-1])
ax.grid()

ax = axes[1, 0]
tdf.plot(ax=ax, style='-k')
ax.set_ylabel('Temperature (degC)')
ax.set_xlim(qdfh.index[0], qdf.index[-1])
ax.grid()

plt.show()