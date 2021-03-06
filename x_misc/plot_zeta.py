"""
Plot results of get_zeta.py

"""

import pandas as pd
import matplotlib.pyplot as plt

import os
import sys
pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import zfun
Ldir = Lfun.Lstart()

# load file
indir = Ldir['LOo'] + 'misc/'
fn = indir + 'zeta_df.p'
zdf = pd.read_pickle(fn)

zdf['sog_minus_jdf'] = zdf['z_sog'] - zdf['z_jdf']
zdf['sog_minus_off'] = zdf['z_sog'] - zdf['z_off']

# plot results
plt.close('all')
zdf.plot(grid=True)
plt.show()
