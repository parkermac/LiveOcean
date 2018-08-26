"""
Code to plot a map for planning TEF sections.
"""
# setup
import matplotlib.pyplot as plt

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

plt.close('all')

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
pfun.add_coast(ax)
pfun.dar(ax)
ax.set_xlim(-126, -122)
ax.set_ylim(47, 50.3)
ax.grid(True)

plt.show()