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

import tef_fun
from importlib import reload
reload(tef_fun)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

plt.close('all')
fig = plt.figure(figsize=(12,8))

def plotit(ax, aa, sect_df):
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.grid(True)
    for sn in sect_df.index:
        x0 = sect_df.loc[sn,'x0']
        x1 = sect_df.loc[sn,'x1']
        y0 = sect_df.loc[sn,'y0']
        y1 = sect_df.loc[sn,'y1']
        xx = (x0+x1)/2
        yy = (y0+y1)/2
        if (xx > aa[0]) & (xx < aa[1]) & (yy > aa[2]) & (yy < aa[3]):
            ax.plot([x0,x1], [y0,y1], '-b', linewidth=3, alpha=.5)
            ax.text(xx, yy, sn, fontsize=12, color='r', fontweight='bold')
    

ax = fig.add_subplot(121)
#aa = [-125.5, -122, 46, 50.3]
aa = [-127.1, -122, 44.5, 50.3]
plotit(ax, aa, sect_df)

ax = fig.add_subplot(222)
aa = [-123.3, -122, 47, 48.4]
plotit(ax, aa, sect_df)

ax = fig.add_subplot(224)
aa = [-123, -122.5, 47.05, 47.35]
plotit(ax, aa, sect_df)

plt.show()