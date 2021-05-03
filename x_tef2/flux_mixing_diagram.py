"""
Use info from TEF sections to calculate values for the
Geyer and MacCready 2014 parameter space diagram.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
from datetime import datetime, timedelta
import netCDF4 as nc

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
x_indir = indir0 + 'extractions/'

outdir = indir00 + 'sill_dyn_plots/'
Lfun.make_dir(outdir)

# get section definitions
sect_df = tef_fun.get_sect_df()

q_df = pd.DataFrame(index=sect_df.index,
    columns=['Ut','H', 'Ur', 'M', 'Fr', 'Ue', 'DS', 'c', 'Ue_non', 'DS_non'])

dp_tef_df, dp_in_sign = flux_fun.get_fluxes(indir0, 'dp')
dp_Qr = (dp_tef_df['Qin'] + dp_tef_df['Qout']).mean()

# constants
Cd = 2.5e-3
om = 1.4e-4
beta = 7.7e-4
g = 9.8
Socn=32

for sect_name in q_df.index:
    
    tef_df, in_sign = flux_fun.get_fluxes(indir0, sect_name)
    qe = (tef_df['Qin'] - tef_df['Qout']).mean()/2
    DS = (tef_df['Sin'] - tef_df['Sout']).mean()/2
        
    ds = nc.Dataset(x_indir + sect_name + '.nc')
    
    if False:
        H = ds['h'][:].mean() # average depth [m]
    else:
        H = ds['h'][:].max() # max depth [m]
        
    A = ds['DA0'][:].sum() # cross-sectional area [m2]
    
    N0 = np.sqrt(beta*g*Socn/H)
    c = np.sqrt(beta*g*Socn*H)
    
    Ut = (np.pi/2) * tef_df['Qtide'].mean() / A
    Qr = (tef_df['Qin'] + tef_df['Qout']).mean()
    if 'ai' in sect_name:
        # This is meant to "correct" the riverflow driving the exchange
        # at Admiralty Inlet.  Perhaps to do this right I should revisit
        # the Sutherland calculation of the fraction of Skagit water that
        # goes through Deception Pass.
        Qr += dp_Qr
    Ur = np.abs(Qr) / A
    
    q_df.loc[sect_name,'H'] = H
    q_df.loc[sect_name,'Ut'] = Ut
    q_df.loc[sect_name,'Ur'] = Ur
    
    Ue = np.abs(qe) / A
    DS = np.abs(DS)
    q_df.loc[sect_name,'Ue'] = Ue
    q_df.loc[sect_name,'DS'] = DS
    
    # derived quantities
    M2 = (Cd*Ut*Ut)/(om*N0*H*H)
    M = np.sqrt(M2)
    
    Fr = Ur / c
    
    q_df.loc[sect_name,'M'] = M
    q_df.loc[sect_name,'Fr'] = Fr
    q_df.loc[sect_name,'c'] = c
    
    # non-simansional versions
    q_df.loc[sect_name,'Ue_non'] = Ue / c
    q_df.loc[sect_name,'DS_non'] = DS / Socn
    
    
# plotting
plt.close('all')
fs = 16
plt.rc('font', size=fs)

# create the GM14 mixing diagram

fig = plt.figure(figsize=(11,11))
ax = fig.add_subplot(111)

M_vec = q_df['M'].to_numpy(dtype=float)
Fr_vec = q_df['Fr'].to_numpy(dtype=float)
DS_vec = q_df['DS'].to_numpy(dtype=float)
Ue_vec = q_df['Ue'].to_numpy(dtype=float)
c_vec = q_df['c'].to_numpy(dtype=float)

ax.scatter(M_vec, Fr_vec , s=100*DS_vec, c=np.log(Ue_vec/c_vec), marker='o', cmap='jet', alpha=.5)
ax.set_yscale('log')
ax.set_xscale('log')

for sn in q_df.index:
    ax.text(q_df.loc[sn,'M'], q_df.loc[sn, 'Fr'], sn, size=fs/2)

# add a line
ax.plot([.3, 1.5], [1e-4,1], '-g', lw=2)

ax.grid(True)
ax.set_xlabel('M')
ax.set_ylabel('Fr')

# plot result vs M and Fr

# remove deception pass
q1_df = q_df.copy()
q1_df.drop('dp', inplace=True)

M_vec = q1_df['M'].to_numpy(dtype=float)
Fr_vec = q1_df['Fr'].to_numpy(dtype=float)
DS_non_vec = q1_df['DS_non'].to_numpy(dtype=float)
Ue_non_vec = q1_df['Ue_non'].to_numpy(dtype=float)

fig = plt.figure(figsize=(11,11))


ax = fig.add_subplot(221)
ax.plot(Fr_vec, DS_non_vec, 'ob', alpha=.5)
ax.set_xlabel('Fr')
ax.set_ylabel('DeltaS/Socn')

ax = fig.add_subplot(222)
ax.plot(Fr_vec, Ue_non_vec, 'ob', alpha=.5)
ax.set_xlabel('Fr')
ax.set_ylabel('Ue/c')

ax = fig.add_subplot(223)
ax.plot(M_vec, DS_non_vec, 'ob', alpha=.5)
ax.set_xlabel('M')
ax.set_ylabel('DeltaS/Socn')

ax = fig.add_subplot(224)
ax.plot(M_vec, Ue_non_vec, 'ob', alpha=.5)
ax.set_xlabel('M')
ax.set_ylabel('Ue/c')


plt.show()
plt.rcdefaults()
    
    
    

