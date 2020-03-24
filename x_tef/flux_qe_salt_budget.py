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

# get section definitions (for dx in dSbar_dx)
sect_df = tef_fun.get_sect_df()

# Note that we trim off the first list item because we are doing our
# dynamical calculations at something like the middle of a sill,
# e.g. ai2 where we can sensible calculate gradients like dSbar_dx.
#
# There will be some dynamical inconsistency for Puget Sound
# because we include Deception Pass.  Similar issue for Salish Sea
# because of Johnstone Strait.
if which_vol == 'Salish Sea':
    seg_list = list(v_lp_df.columns)
    seg_list = seg_list[1:]
    sect_sign_dict = {'jdf2':1, 'sog5':-1}
elif which_vol == 'Puget Sound':
    seg_list = (flux_fun.ssA + flux_fun.ssM + flux_fun.ssT
        + flux_fun.ssS + flux_fun.ssW + flux_fun.ssH)
    seg_list = seg_list[1:]
    sect_sign_dict = {'ai2':1, 'dp':1}
elif which_vol == 'Hood Canal':
    seg_list = flux_fun.ssH
    seg_list = seg_list[1:]
    sect_sign_dict = {'hc2':1}
elif which_vol == 'South Sound':
    seg_list = flux_fun.ssT + flux_fun.ssS
    seg_list = seg_list[1:]
    sect_sign_dict = {'tn2':1}

sv_lp_df = sv_lp_df[seg_list]
sv_lp_df = sv_lp_df[seg_list]

river_list = []
for seg_name in seg_list:
    seg = flux_fun.segs[seg_name]
    river_list = river_list + seg['R']
    riv_df = pd.read_pickle(Ldir['LOo'] + 'river/'
        + Ldir['gtag'] + '_'+year_str+'.01.01_'+year_str+'.12.31.p')
riv_df.index += timedelta(days=0.5)
riv_df = riv_df.loc[sv_lp_df.index, river_list]
    
def get_fluxes(sect_name, in_sign=1):
    # form time series of net 2-layer transports into (+) and out of (-) the volume
    bulk = pickle.load(open(indir0 + 'bulk/' + sect_name + '.p', 'rb'))
    QQ = bulk['QQ']
    SS = bulk['SS']
    ot = bulk['ot']
    dt2 = []
    for tt in ot:
        dt2.append(Lfun.modtime_to_datetime(tt))
    # separate inflowing and outflowing transports
    QQp = QQ.copy()
    QQp[QQ<=0] = np.nan
    QQm = QQ.copy()
    QQm[QQ>=0] = np.nan
    # form two-layer versions of Q and S
    if in_sign == 1:
        Qin = np.nansum(QQp, axis=1)
        QSin = np.nansum(QQp*SS, axis=1)
        Qout = np.nansum(QQm, axis=1)
        QSout = np.nansum(QQm*SS, axis=1)
    elif in_sign == -1:
        Qin = -np.nansum(QQm, axis=1)
        QSin = -np.nansum(QQm*SS, axis=1)
        Qout = -np.nansum(QQp, axis=1)
        QSout = -np.nansum(QQp*SS, axis=1)
    Sin = QSin/Qin
    Sout = QSout/Qout
    fnet = bulk['fnet_lp'] # net tidal energy flux
    qabs = bulk['qabs_lp']
    tef_df = pd.DataFrame(index=dt2)#, columns=['Qin','Qout','QSin','QSout','Sin','Sout','fnet'])
    tef_df.loc[:,'Qin']=Qin
    tef_df.loc[:,'Qout']=Qout
    tef_df.loc[:,'QSin']=QSin
    tef_df.loc[:,'QSout']=QSout
    tef_df.loc[:,'Sin']=Sin
    tef_df.loc[:,'Sout']=Sout
    tef_df.loc[:,'fnet'] = fnet * in_sign
    tef_df.loc[:,'Qabs'] = qabs

    return tef_df
    
tef_df_dict = {}
for sn in sect_sign_dict.keys():
    in_sign = sect_sign_dict[sn]
    tef_df_dict[sn] = get_fluxes(sn, in_sign=in_sign)
    
# hack to get dsdx
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
    
df1 = get_fluxes(sea_sect)
df3 = get_fluxes(land_sect)

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

# various things for the dynamicsl scalings
dSbar_dx = ((df1['Sin']+df1['Sout'])/2-(df3['Sin']+df3['Sout'])/2)/DX
DF = df1['fnet']-df3['fnet'] # Net loss of tidal energy flux in region
A0 = v_df.loc[sea_seg,'area m2'] + v_df.loc[land_seg,'area m2']
V0 = v_df.loc[sea_seg,'volume m3'] + v_df.loc[land_seg,'volume m3']
H0 = V0/A0 # average depth of surrounding segments
B0 = V0/(H0*DX)

# volume budget
vol_df = pd.DataFrame(index=sv_lp_df.index, columns=['Qin','-Qout', 'Qabs', 'Qr','dV_dt','dV_dt_check'])
vol_df.loc[:,['Qin','-Qout','Qabs']] = 0
for sect_name in sect_sign_dict.keys():
    df = tef_df_dict[sect_name]
    vol_df.loc[:,'Qin'] = vol_df.loc[:,'Qin'] + df.loc[:,'Qin']
    vol_df.loc[:,'-Qout'] = vol_df.loc[:,'-Qout'] - df.loc[:,'Qout']
    vol_df.loc[:,'Qabs'] = vol_df.loc[:,'Qabs'] + df.loc[:,'Qabs']
qr = riv_df.sum(axis=1)
vol_df.loc[:,'Qr'] = qr
v = v_lp_df[seg_list].sum(axis=1).values
vol_df.loc[1:-1, 'dV_dt'] = (v[2:] - v[:-2]) / (2*86400)

# salt budget
salt_df = pd.DataFrame(index=sv_lp_df.index, columns=['QSin','-QSout','dSnet_dt','fnet'])
salt_df.loc[:,['QSin','-QSout','fnet']] = 0
for sect_name in sect_sign_dict.keys():
    df = tef_df_dict[sect_name]
    salt_df.loc[:,'QSin'] = salt_df.loc[:,'QSin'] + df.loc[:,'QSin']
    salt_df.loc[:,'-QSout'] = salt_df.loc[:,'-QSout'] - df.loc[:,'QSout']
    salt_df.loc[:,'fnet'] = salt_df.loc[:,'fnet'] + df.loc[:,'fnet']
sn = sv_lp_df[seg_list].sum(axis=1).values
salt_df.loc[1:-1, 'dSnet_dt'] = (sn[2:] - sn[:-2]) / (2*86400)
# add a few more columns to plot in a different way
salt_df['Qe'] = (vol_df['Qin'] + vol_df['-Qout'])/2
salt_df['Qr'] = (vol_df['-Qout'] - vol_df['Qin'])
salt_df['Sin'] = salt_df['QSin']/vol_df['Qin']
salt_df['Sout'] = salt_df['-QSout']/vol_df['-Qout']
salt_df['DS'] = salt_df['Sin'] - salt_df['Sout']
salt_df['Sbar'] = (salt_df['QSin']/vol_df['Qin'] + salt_df['-QSout']/vol_df['-Qout'])/2
salt_df['QeDS'] = salt_df['Qe'] * salt_df['DS']
salt_df['-QrSbar'] = -salt_df['Qr'] * salt_df['Sbar']

# salt budget in exchange-flow terms
sf_df = salt_df[['QeDS','-QrSbar','dSnet_dt']]
sf_df['Error'] = sf_df['dSnet_dt'] - sf_df['QeDS'] - sf_df['-QrSbar']
# convert to 1000 g/kg m3/s
sf_df = sf_df/1000

# make sure everything is numeric
for cn in vol_df.columns:
    vol_df[cn] = pd.to_numeric(vol_df[cn])
for cn in salt_df.columns:
    salt_df[cn] = pd.to_numeric(salt_df[cn])
for cn in sf_df.columns:
    sf_df[cn] = pd.to_numeric(sf_df[cn])

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
dyn_df['Qe_pred'] = (Ue * B0 * H0 / 4)/1000
dyn_df['Qe'] = salt_df['Qe']/1000
dyn_df['DS_pred'] = H0**2 * dSbar_dx * Ue / (12 * Ks)
dyn_df['DS'] = salt_df['DS']
dyn_df['QeDS_pred'] = dyn_df['Qe_pred'] * dyn_df['DS_pred']
dyn_df['QeDS'] = sf_df['QeDS']
dyn_df['Qabs'] = (vol_df['Qabs']/1000)/7

for cn in dyn_df.columns:
    dyn_df[cn] = pd.to_numeric(dyn_df[cn])

plt.close('all')
fig = plt.figure(figsize=(20,10))

ax = fig.add_subplot(211)
sf_df[['QeDS','-QrSbar','dSnet_dt','Error']].plot(ax=ax, grid=True, legend=False)
ax.set_title('%s Salt Budget $(10^{3}\ g\ kg^{-1}\ m^{3}s^{-1})$' % (which_vol))
ax.set_xlim(sf_df.index[0], sf_df.index[-1])
legh = ax.legend(labels=['$Q_e \Delta S$', '$-Q_R S_{bar}$', '$d S_{net} / dt$','Error'])

# # get some things for correlations
# QeDS = pd.to_numeric(sf_df['QeDS']).to_numpy()
# dSdx = pd.to_numeric(sf_df['dSdx']).to_numpy()
# Dissipation = pd.to_numeric(sf_df['Dissipation']).to_numpy()
# Forcing = pd.to_numeric(sf_df['Forcing']).to_numpy()
# rf = np.corrcoef(QeDS, Forcing)[0,1]
# rd = np.corrcoef(QeDS, Dissipation)[0,1]
# rs = np.corrcoef(QeDS, dSdx)[0,1]
# rdict = {'Forcing':rf, 'Dissipation':rd, 'dSdx':rs}

ax = fig.add_subplot(245)
dyn_df.plot(x='Qe', y='Qe_pred', ax=ax, grid=True, style='og', markersize=2)
ax.axis('square')
ax.axis([0,60,0,60])

ax = fig.add_subplot(246)
dyn_df.plot(x='DS', y='DS_pred', ax=ax, grid=True, style='og', markersize=2)
ax.axis('square')
ax.axis([0,3,0,3])

ax = fig.add_subplot(247)
dyn_df.plot(x='QeDS', y='QeDS_pred', ax=ax, grid=True, style='og', markersize=2)
ax.axis('square')
ax.axis([0,100,0,100])


ax = fig.add_subplot(248)
dyn_df.plot(x='Qe', y='Qabs', ax=ax, grid=True, style='og', markersize=2)
ax.axis('square')
ax.axis([0,100,0,100])


plt.show()
