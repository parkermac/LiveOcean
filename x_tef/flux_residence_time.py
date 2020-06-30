"""
Designed to calculate residence time from all of the "IC" experiments.

Also makes a plot.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
Ldir['gtagex'] = Ldir['gtag'] + '_lo8b'
import zfun

import tef_fun
import flux_fun

# Input directory
indir0 = Ldir['LOo'] + 'tef/'
outdir = indir0 + 'misc_figs_cas6/'
indir = indir0 + 'flux_engine/' + Ldir['gtagex'] + '/'
voldir = indir0 + 'volumes_' + Ldir['gridname'] + '/'

# load a Series of the volumes of each segment, created by flux_get_vol.py
v_df = pd.read_pickle(voldir + 'volumes.p')
V = flux_fun.get_V(v_df)

ic_list = [item for item in os.listdir(indir) if item[:3]=='IC_']
ic_list.sort()

basin_list = ['Salish', 'SoG', 'PS', 'Whidbey', 'HoodCanal', 'HoodCanalInner', 'SouthSound']
basin_str = ['Salish Sea', 'Strait of Georgia', 'Puget Sound',
    'Whidbey Basin', 'Hood Canal', 'Inner Hood Canal', 'South Sound']
basin_dict = dict(zip(basin_list,basin_str))

# get list of experiments
exp_list = []
for infile in ic_list:
    # infile is like "'IC_Whidbey_2019_winter.p'"
    exp = infile.replace('IC_','')
    exp = exp.replace('.p','')
    exp_list.append(exp)
    
# add an x value for each experiment
x = np.arange(1,len(basin_list)+1,dtype='float') - .5
x_dict = dict(zip(basin_list,x))

# and an offset to x based on year
year_list = ['2017','2018','2019']
year_xlist = [-.25, 0, .25]
year_xdict = dict(zip(year_list, year_xlist))

tres_df = pd.DataFrame(0, index = exp_list,columns=['basin','year', 'season', 'tres', 'tres0','x'])
        # tres0 is the residence time without reflux


for infile in ic_list:
    # infile is like "'IC_Whidbey_2019_winter.p'"
    
    source = infile.split('_')[0] + '_' + infile.split('_')[1]
    
    
    
    # get itemized information for the DataFrame
    exp = infile.replace('IC_','')
    exp = exp.replace('.p','')
    basin, year, season = exp.split('_')
    
    seg2_list = flux_fun.ic_seg2_dict[source]
        
    aa = pd.read_pickle(indir + infile)
    
    this_aa = aa.loc[:,seg2_list]
    this_V = V[seg2_list]
    net_V = this_V.sum()

    this_net_aa = this_aa.copy()

    for sn in this_V.index:
        VV = this_V[sn]
        this_net_aa.loc[:,sn] = this_net_aa.loc[:,sn] * VV
    
    # make a Series of mean concentration in the volume
    mean_c = this_net_aa.sum(axis=1) / net_V
    
    # find e-folding time
    td = mean_c.index.to_numpy()
    mc = mean_c.to_numpy()
    ind_ef = np.argwhere(mc < 1/np.e)[0]
    tres = td[ind_ef] # residence time in days
    
    # also load the A matrix to allow us to calculate the "unrefluxed"
    # residence time
    q_df = pd.read_pickle(Ldir['LOo'] + 'tef/' +
        Ldir['gtagex'] + '_' + year + '.01.01_' + year + '.12.31/' +
        'flux/q_df_' + season + '.p')
    if basin == 'HoodCanalInner':
        qin = q_df.loc['H3_s','H2_s']
    elif basin == 'HoodCanal':
        qin = q_df.loc['H1_s','A3_s']
    elif basin == 'SouthSound':
        qin = q_df.loc['S1_s','T2_s']
    elif basin == 'PS':
        qin = q_df.loc['A1_s','J4_s'] + q_df.loc['W4_s','J4_s']
    elif basin == 'Whidbey':
        qin = q_df.loc['W1_s','M1_s'] + q_df.loc['W4_s','J4_s']
    elif basin == 'Salish':
        qin = q_df.loc['J1_s','ocean_s'] + q_df.loc['G6_s','ocean_s']
    elif basin == 'SoG':
        qin = q_df.loc['G1_s','J4_s'] + q_df.loc['G6_s','ocean_s']
    else:
        print(basin)
        print('unsupported qin')
        print('')
    tres0 = (net_V/qin)/86400
    
    # store results
    
    xx = x_dict[basin] + year_xdict[year]
    tres_df.loc[exp,['basin', 'year','season','tres','tres0','x']] = [basin, year, season, tres, tres0,xx]
            
    if False:
        print(' %10s: Tres = %0.2f days' % (exp, tres))


# residence time figure
plt.close('all')
fs = 16
plt.rc('font', size=fs)

fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)

season_list = ['winter','spring','summer','fall']
#season_clist = ['cornflowerblue', 'limegreen', 'lightsalmon', 'goldenrod']
season_clist = ['b', 'g', 'r', 'orange']
season_cdict = dict(zip(season_list, season_clist))

for season in season_list:
    df = tres_df[tres_df['season']==season]
    df.plot(x='x', y='tres', ax=ax, marker='o',
        color=season_cdict[season], markersize=24, alpha=.7,
        linestyle='None', legend=False)

for season in season_list:
    df = tres_df[tres_df['season']==season]
    df.plot(x='x', y='tres0', ax=ax, marker='*',
        color=season_cdict[season], markersize=12, alpha=.7,
        linestyle='None', legend=False)
        
ax.set_xlim(0,7)
ax.set_xticks(x)
ax.set_xlabel('')

ax.set_ylabel('Residence Time [days]')
ax.set_ylim(0,600)

for X in np.arange(1,7):
    ax.axvline(x=X, color='gray')
for Y in np.arange(1,6)*100:
    ax.axhline(y=Y, color='gray', linestyle='--')

ax.set_xticklabels(basin_str, rotation=70, weight='bold', size=1.5*fs)

ax.text(2.25, 250, '2017', rotation=90, style='italic', ha='center', va='center', weight='bold', size=fs*1.3)
ax.text(2.5, 250, '2018', rotation=90, style='italic', ha='center', va='center', weight='bold', size=fs*1.3)
ax.text(2.75, 250, '2019', rotation=90, style='italic', ha='center', va='center', weight='bold', size=fs*1.3)

ii = 0
for season in season_list:
    ax.text(2.1, 575 - 50*ii, season.title(),
        color=season_cdict[season], weight='bold', size=fs*1.3, va='center')
    ii += 1

fig.tight_layout()
fig.savefig(outdir + 'all_residence_times.png')
plt.show()
plt.rcdefaults()

