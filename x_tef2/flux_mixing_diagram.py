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

testing = False
if testing:
    sect_list = ['ai1']
else:
    sect_list = list(sect_df.index)

# initialize DataFrame
q_df = pd.DataFrame(index=sect_list,
    columns=['Ut', 'Qprism','H', 'Qr','Ur', 'M', 'Fr', 'Qe', 'Ue', 'DS', 'c', 'Ue_non', 'DS_non'])

# get Socn
jdf1_tef_df, jdf1_in_sign = flux_fun.get_fluxes(indir0, 'jdf1')
Socn = jdf1_tef_df['Sin'].max()
print('Socn = %0.2f' % (Socn))

# constants
Cd = 2.5e-3
om = 1.4e-4
beta = 7.7e-4
g = 9.8

for sect_name in sect_list:
    
    tef_df, in_sign = flux_fun.get_fluxes(indir0, sect_name)
    Qe = (tef_df['Qin'] - tef_df['Qout']).mean()/2
    DS = (tef_df['Sin'] - tef_df['Sout']).mean()
        
    ds = nc.Dataset(x_indir + sect_name + '.nc')
    H = ds['h'][:].max() # max depth [m]
    A = ds['DA0'][:].sum() # cross-sectional area [m2]
    
    N0 = np.sqrt(beta*g*Socn/H)
    c = np.sqrt(beta*g*Socn*H)
    
    Qprism = tef_df['Qtide'].mean() / 2
    Ut = (np.pi/2) * tef_df['Qtide'].mean() / A
    
    # use Freshwater Flux as an alternate way to calculate Qr and Ur
    Qr = -( tef_df['Qin']*(Socn-tef_df['Sin']) + tef_df['Qout']*(Socn-tef_df['Sout']) ).mean()/Socn
    
    do_sect = True
    if Qr < 1:
        print('Dropping Section: Qr negative for ' + sect_name)
        do_sect = False
        q_df = q_df.drop(sect_name)
        
    if do_sect:
        
        q_df.loc[sect_name,'Qprism'] = Qprism
        
        Ur = np.abs(Qr/A)
    
        if Qe < 0:
            print('Qe negative for ' + sect_name)
        Ue = np.abs(Qe/A) # should we use A/2?
    
        q_df.loc[sect_name,'H'] = H
        q_df.loc[sect_name,'Ut'] = Ut
        q_df.loc[sect_name,'Qe'] = Qe
        q_df.loc[sect_name,'Qr'] = Qr
        q_df.loc[sect_name,'Ur'] = Ur
    
        q_df.loc[sect_name,'Ue'] = Ue
        q_df.loc[sect_name,'DS'] = DS
    
        # derived quantities
        M2 = (Cd*Ut*Ut)/(om*N0*H*H)
        M = np.sqrt(M2)
    
        Fr = Ur / c
    
        q_df.loc[sect_name,'M'] = M
        q_df.loc[sect_name,'Fr'] = Fr
        q_df.loc[sect_name,'c'] = c
    
        # non-dimensional versions
        q_df.loc[sect_name,'Ue_non'] = Ue / c
        q_df.loc[sect_name,'DS_non'] = DS / Socn

plt.close('all')
    
if True:
    
    # plotting
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

    fig.savefig(outdir + 'MG14_diagram.png')
    plt.show()
    plt.rcdefaults()

if True:
    # plot result vs M and Fr

    # plotting
    fs = 12
    plt.rc('font', size=fs)

    if True:
        # remove deception pass
        print('Dropping Deception Pass')
        print(q_df.loc['dp',:])
        q_df = q_df.drop('dp')

    fig = plt.figure(figsize=(14,11))

    ax = fig.add_subplot(221)
    q_df.plot(x='Fr', y ='DS_non', style='ob', alpha=.5, legend=False, ax=ax)
    ax.set_xlabel('Fr')
    ax.set_ylabel('DeltaS/Socn')

    ax = fig.add_subplot(222)
    q_df.plot(x='Fr', y ='Ue_non', style='ob', alpha=.5, legend=False, ax=ax)
    ax.set_xlabel('Fr')
    ax.set_ylabel('Ue/c')
    for sn in q_df.index:
        ax.text(q_df.loc[sn,'Fr'], q_df.loc[sn, 'Ue_non'], sn, size=fs/2)

    ax = fig.add_subplot(223)
    q_df.plot(x='M', y ='DS_non', style='ob', alpha=.5, legend=False, ax=ax)
    ax.set_xlabel('M')
    ax.set_ylabel('DeltaS/Socn')

    ax = fig.add_subplot(224)
    q_df.plot(x='M', y ='Ue_non', style='ob', alpha=.5, legend=False, ax=ax)
    ax.set_xlabel('M')
    ax.set_ylabel('Ue/c')
    
    fig.savefig(outdir + 'MG14_property_property.png')
    plt.show()
    plt.rcdefaults()
    
if True:
    # plot exchange vs. tide in various ways

    # plotting
    fs = 14
    plt.rc('font', size=fs)

    if True and 'dp' in q_df.index:
        # remove deception pass
        print('Dropping Deception Pass')
        print(q_df.loc['dp',:])
        q_df = q_df.drop('dp')

    fig = plt.figure(figsize=(14,11))

    ax = fig.add_subplot(221)
    q_df.plot(x='Qprism', y ='Qe', style='o', alpha=.5, legend=False, ax=ax, loglog=True)
    for sn in q_df.index:
        ax.text(q_df.loc[sn,'Qprism'], q_df.loc[sn, 'Qe'], sn, size=fs/2)
    ax.set_xlabel('Qprism')
    ax.set_ylabel('Qe')

    ax = fig.add_subplot(222)
    q_df.plot(x='Qprism', y ='Qe', style='o', alpha=.5, legend=False, ax=ax)
    ax.set_xlabel('Qprism')
    ax.set_ylabel('Qe')

    ax = fig.add_subplot(223)
    q_df.plot(x='Ut', y ='Ue', style='o', alpha=.5, legend=False, ax=ax, loglog=True)
    for sn in q_df.index:
        ax.text(q_df.loc[sn,'Ut'], q_df.loc[sn, 'Ue'], sn, size=fs/2)
    ax.set_xlabel('Ut')
    ax.set_ylabel('Ue')

    ax = fig.add_subplot(224)
    q_df.plot(x='Ut', y ='Ue', style='o', alpha=.5, legend=False, ax=ax)
    ax.set_xlabel('Ut')
    ax.set_ylabel('Ue')
    
    fig.savefig(outdir + 'MG14_exchange_vs_tide.png')
    plt.show()
    plt.rcdefaults()    
    

