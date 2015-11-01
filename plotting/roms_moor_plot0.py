"""
Plots mooring records.
"""

# setup
import os; import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path: sys.path.append(alp)
import Lfun; reload(Lfun)
import numpy as np
import zfun; reload(zfun) # plotting functions
import matfun; reload(matfun) # functions for working with mat files

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'moor/'

# choose the type of plot to make
print '\n%s\n' % '** Choose mooring file to plot **'
m_list_raw = os.listdir(indir)
m_list = []
for m in m_list_raw:
    if m[-2:] == '.p':
        m_list.append(m)
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print str(npt) + ': ' + m_list[npt]
my_npt = int(raw_input('-- Input number -- '))
inname = m_dict[my_npt]
  
import cPickle as pickle
V, v1_list, v2_list, v3_list, G = pickle.load( open( indir + inname, 'rb' ) )

import matplotlib.pyplot as plt
plt.close()
NR = 3; NC = 4
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,8), squeeze=False)

days = (V['ocean_time'] - V['ocean_time'][0])/86400.
cc = 0
for vn in v3_list:
    ir = int(np.floor(cc/NC))
    ic = int(cc - NC*ir)
    ax = axes[ir, ic]
    ax.plot(days, V[vn][-1,:])
    ax.set_xlim(0, days.max())
    ax.ticklabel_format(useOffset=False, axis='y')
    ax.set_title(vn)
    cc += 1
for vn in v2_list:
    ir = int(np.floor(cc/NC))
    ic = int(cc - NC*ir)
    ax = axes[ir, ic]
    ax.plot(days, V[vn])
    ax.set_xlim(0, days.max())
    ax.ticklabel_format(useOffset=False, axis='y')
    ax.set_title(vn)
    cc += 1
plt.show()
