"""
Compare exchange flow and tidal transport for ALL
sections.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
from datetime import datetime, timedelta

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zrfun
import zfun

import tef_fun
import flux_fun

from time import time

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
# associated with lines like QQp[QQ<=0] = np.nan

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', type=str, default='cas6')
parser.add_argument('-t', '--tag', type=str, default='v3')
parser.add_argument('-x', '--ex_name', type=str, default='lo8b')
parser.add_argument('-y', '--year', type=int, default=2017)
args = parser.parse_args()
year_str = str(args.year)

# Get Ldir
Ldir = Lfun.Lstart(args.gridname, args.tag)
gtagex = args.gridname + '_' + args.tag + '_' + args.ex_name

# select input/output location
run_name = gtagex+'_'+year_str+'.01.01_'+year_str+'.12.31'
indir00 = Ldir['LOo'] + 'tef2/'
indir0 = indir00 + run_name + '/'
indir = indir0 + 'flux/'

outdir = indir00 + 'sill_dyn_plots/'
Lfun.make_dir(outdir)

# get section definitions
sect_df = tef_fun.get_sect_df()


q_df = pd.DataFrame(index=sect_df.index, columns=['Qe','Qprism','Ftide'])

for sect_name in q_df.index:
    tef_df, in_sign = flux_fun.get_fluxes(indir0, sect_name)
    qe = (tef_df['Qin'] - tef_df['Qout'])/2
    QE = qe.mean()/1000 # convert to 1000 m3/s
    QT = tef_df['Qtide'].mean()/1000 # convert to 1000 m3/s
    FT = tef_df['Ftide'].mean() # units are MW
    q_df.loc[sect_name,'Qe'] = QE
    q_df.loc[sect_name,'Qprism'] = QT/2
    q_df.loc[sect_name,'Ftide'] = FT

# fit a line
x = q_df['Qprism'].to_numpy(dtype=float)
y = q_df['Qe'].to_numpy(dtype=float)
BB = np.polyfit(x,y,1)
dydx = BB[0] # slope
y0 = 0 #BB[1] # y-intercept [force to zero]

xx = np.linspace(.1,x.max(),1000)
yy = y0 + dydx*xx

plt.close('all')

fs=18
plt.rc('font', size=fs)

# linear scale for axes
fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111)
q_df.plot(x='Qprism', y = 'Qe', style='ok', ax=ax, grid=True, legend=False, ms=14, alpha=.2)
ax.set_xlim(xmin=0)
ax.set_ylim(ymin=0)
ax.set_xlabel(r'Annual mean of Qprism $[10^{3}m^{3}s^{-1}]$')
ax.set_ylabel(r'Annual mean of Qe $[10^{3}m^{3}s^{-1}]$')
# add the fit line
ax.plot(xx,yy,'-c', lw=2)
ax.text(.05,.9,'y = %0.2f*x + %0.2f' % (dydx, y0), transform=ax.transAxes, color='c', weight='bold')
# add labels
for sect_name in q_df.index:
    ax.text(q_df.loc[sect_name,'Qprism'],q_df.loc[sect_name,'Qe'],sect_name,
    ha='center', va='center', size=.7*fs, rotation=-45)

# log scale for axes
fig2 = plt.figure(figsize=(12,12))
ax = fig2.add_subplot(111)
q_df.plot(x='Qprism', y = 'Qe', style='ok', ax=ax, loglog=True, grid=True, legend=False, ms=14, alpha=.2)
ax.set_xlim(xmin=.1)
ax.set_ylim(ymin=.1)
ax.set_xlabel(r'Annual mean of Qprism $[10^{3}m^{3}s^{-1}]$')
ax.set_ylabel(r'Annual mean of Qe $[10^{3}m^{3}s^{-1}]$')
# add the fit line
ax.loglog(xx,yy,'-c', lw=2)
ax.text(.05,.9,'y = %0.2f*x + %0.2f' % (dydx, y0), transform=ax.transAxes, color='c', weight='bold')
# add labels
for sect_name in q_df.index:
    ax.text(q_df.loc[sect_name,'Qprism'],q_df.loc[sect_name,'Qe'],sect_name,
    ha='center', va='center', size=.7*fs, rotation=-45)

plt.show()

plt.rcdefaults()

fig.savefig(outdir + 'allmean_' + year_str + '.png')
fig2.savefig(outdir + 'allmean_' + year_str + '_loglog.png')
