"""
Calculates an exchange-flow oriented salt budget from TEF terms, and
explores dynamical scaling:

does Qe behave as expected relative to dSbar_dx and K?

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
from importlib import reload
reload(flux_fun)

from time import time

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
# associated with lines like QQp[QQ<=0] = np.nan

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', type=str, default='cas6')
parser.add_argument('-t', '--tag', type=str, default='v3')
parser.add_argument('-x', '--ex_name', type=str, default='lo8b')
parser.add_argument('-y', '--year', type=int, default=2018)
parser.add_argument('-v', '--volume', type=str, default='Puget Sound')
args = parser.parse_args()
which_vol = args.volume
year_str = str(args.year)

# Get Ldir
Ldir = Lfun.Lstart(args.gridname, args.tag)
gtagex = args.gridname + '_' + args.tag + '_' + args.ex_name

# select input/output location
run_name = gtagex+'_'+year_str+'.01.01_'+year_str+'.12.31'
indir00 = Ldir['LOo'] + 'tef/'
indir0 = indir00 + run_name + '/'
indir = indir0 + 'flux/'

# load low passed segment volume and net salt DataFrames
v_lp_df = pd.read_pickle(indir + 'daily_segment_volume.p')
sv_lp_df = pd.read_pickle(indir + 'daily_segment_net_salt.p')

# get volumes
voldir = indir00 + 'volumes_' + Ldir['gridname'] + '/'
v_df = pd.read_pickle(voldir + 'volumes.p')

# get section definitions
sect_df = tef_fun.get_sect_df()

# Note that we trim off the first list item because we are doing our
# dynamical calculations at something like the middle of a sill,
# e.g. ai2 where we can sensibly calculate gradients like dSbar_dx.
if which_vol == 'Salish Sea':
    seg_list = list(v_lp_df.columns)
    seg_list = seg_list[1:]
    sect_sign_dict = {'jdf2':1}
elif which_vol == 'Puget Sound':
    seg_list = (flux_fun.ssA + flux_fun.ssM + flux_fun.ssT
        + flux_fun.ssS + flux_fun.ssW + flux_fun.ssH)
    seg_list = seg_list[1:]
    sect_sign_dict = {'ai2':1}
elif which_vol == 'Hood Canal':
    seg_list = flux_fun.ssH
    seg_list = seg_list[1:]
    sect_sign_dict = {'hc2':1}
elif which_vol == 'South Sound':
    seg_list = flux_fun.ssT + flux_fun.ssS
    seg_list = seg_list[1:]
    sect_sign_dict = {'tn2':1}
elif which_vol == 'Strait of Georgia':
    seg_list = flux_fun.ssG
    seg_list = seg_list[1:]
    sect_sign_dict = {'sji2':1}

v_lp_df = v_lp_df[seg_list]
sv_lp_df = sv_lp_df[seg_list]

river_list = []
for seg_name in seg_list:
    seg = flux_fun.segs[seg_name]
    river_list = river_list + seg['R']
    riv_df = pd.read_pickle(Ldir['LOo'] + 'river/'
        + Ldir['gtag'] + '_'+year_str+'.01.01_'+year_str+'.12.31.p')
riv_df.index += timedelta(days=0.5)
riv_df = riv_df[river_list]
        
tef_df_dict = {}
for sn in sect_sign_dict.keys():
    in_sign = sect_sign_dict[sn]
    tef_df_dict[sn] = flux_fun.get_fluxes(indir0, sn, in_sign=in_sign)
    
vol_df, salt_df, vol_rel_err, salt_rel_err, salt_rel_err_qe = flux_fun.get_budgets(
    sv_lp_df, v_lp_df, riv_df, tef_df_dict, seg_list)
    
# getting gradients across the seaward section
if which_vol == 'Salish Sea':
    sea_sect = 'jdf1'; land_sect = 'jdf3'
    sea_seg = 'J1'; land_seg = 'J2'
elif which_vol == 'Puget Sound':
    sea_sect = 'ai1'; land_sect = 'ai3'
    sea_seg = 'A1'; land_seg = 'A2'
elif which_vol == 'South Sound':
    sea_sect = 'tn1'; land_sect = 'tn3'
    sea_seg = 'T1'; land_seg = 'T2'
elif which_vol == 'Hood Canal':
    sea_sect = 'hc1'; land_sect = 'hc3'
    sea_seg = 'H1'; land_seg = 'H2'
elif which_vol == 'Strait of Georgia':
    sea_sect = 'sji1'; land_sect = 'sog1'
    sea_seg = 'G1'; land_seg = 'G2'
    
df1 = flux_fun.get_fluxes(indir0, sea_sect)
df3 = flux_fun.get_fluxes(indir0, land_sect)

# get DX for dSbar_dx
sea_lon = (sect_df.loc[sea_sect,'x0'] + sect_df.loc[sea_sect,'x1'])/2
sea_lat = (sect_df.loc[sea_sect,'y0'] + sect_df.loc[sea_sect,'y1'])/2
land_lon = (sect_df.loc[land_sect,'x0'] + sect_df.loc[land_sect,'x1'])/2
land_lat = (sect_df.loc[land_sect,'y0'] + sect_df.loc[land_sect,'y1'])/2
mean_lon = (sea_lon + land_lon)/2
mean_lat = (sea_lat + land_lat)/2
sea_x, sea_y = zfun.ll2xy(sea_lon, sea_lat, mean_lon, mean_lat)
land_x, land_y = zfun.ll2xy(land_lon, land_lat, mean_lon, mean_lat)
DX = np.sqrt((sea_x-land_x)**2 + (sea_y-land_y)**2)

# various things for the dynamical scalings
dSbar_dx = ((df1['Sin']+df1['Sout'])/2-(df3['Sin']+df3['Sout'])/2)/DX
DF = df1['Ftide']-df3['Ftide'] # Net loss of tidal energy flux in region
A0 = v_df.loc[sea_seg,'area m2'] + v_df.loc[land_seg,'area m2']
V0 = v_df.loc[sea_seg,'volume m3'] + v_df.loc[land_seg,'volume m3']
H0 = V0/A0 # average depth of surrounding segments
B0 = V0/(H0*DX)

# dynamical scalings
a = 2.5 * 0.028 # the 2.5 is a fudge factor to get Qe to match Qe_pred
Cd = 2.6e-3
h = H0/2 # m
rho = 1027 # kg m-3
g = 9.8 # m2 s-1
beta = 7.7e-4
dyn_df = pd.DataFrame(index=salt_df.index)
Km = (a**3 * Cd**2 * h**3 * DF / (rho*A0))**(1/3)
Ks = Km/2.2
dyn_df['Km'] = Km
dyn_df['Ks'] = Ks
dyn_df['dSbar_dx'] = dSbar_dx
Ue = g * beta * dSbar_dx * H0**3 / (48 * Km)
dyn_df['Ue'] = Ue
dyn_df['Qe_pred'] = (Ue * B0 * H0 / 4)
dyn_df['Qe'] = salt_df['Qe']
dyn_df['DS_pred'] = H0**2 * dSbar_dx * Ue / (12 * Ks)
dyn_df['DS'] = salt_df['DS']
dyn_df['QeDS_pred'] = dyn_df['Qe_pred'] * dyn_df['DS_pred']
dyn_df['QeDS'] = salt_df['QeDS']
dyn_df['Qe_tide_pred'] = (vol_df['Qtide'])/7

for cn in dyn_df.columns:
    dyn_df[cn] = pd.to_numeric(dyn_df[cn])
    if 'Q' in cn:
        # convert transports to 1000 m3/s
        dyn_df[cn] = dyn_df[cn]/1000 

plt.close('all')
fig = plt.figure(figsize=(10,12))

axcount = 1
for vtup in [('Qe','Qe_pred'),('DS','DS_pred'),('QeDS','QeDS_pred'),('Qe','Qe_tide_pred')]:
    ax = fig.add_subplot(2,2,axcount)
    v0 = dyn_df[vtup[0]].to_numpy()
    v1 = dyn_df[vtup[1]].to_numpy()
    r = np.corrcoef(v0, v1)[0,1]
    ax.plot(v0, v1, 'og', markersize=2)
    ax.grid(True)
    ax.axis('square')
    axmax = v0.mean() + 2*v0.std()
    ax.set_xlim(0,axmax)
    ax.set_ylim(0,axmax)
    ax.plot([0,axmax],[0,axmax],'-k')
    ax.set_xlabel(vtup[0])
    ax.set_ylabel(vtup[1])
    ax.text(.95, .1, '$r^{2}=%0.2f$' % (r*r), ha='right', transform = ax.transAxes)
    axcount += 1
fig.suptitle(which_vol)

plt.show()
