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

vftag = '' # use to switch between different volume partition fractions
# '_10' = 10/90
# '' = 20/80 original
# '_30' = 30/70

# Input directory
indir0 = Ldir['LOo'] + 'tef/'
outdir = indir0 + 'misc_figs_cas6/'
indir = indir0 + 'flux_engine/' + Ldir['gtagex'] + vftag + '/'
voldir = indir0 + 'volumes_' + Ldir['gridname'] + '/'

# load a Series of the volumes of each segment, created by flux_get_vol.py
v_df = pd.read_pickle(voldir + 'volumes.p')
V = flux_fun.get_V(v_df)

ic_list = [item for item in os.listdir(indir) if (item[:3]=='IC_') and ('0.p' not in item) and ('full' not in item)]
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

tres_df = pd.DataFrame(0, index = exp_list,columns=['basin','year', 'season', 'tres', 'tres0','tres00','x'])
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
    aa0 = pd.read_pickle(indir + infile.replace('.p','0.p'))
    
    this_aa = aa.loc[:,seg2_list]
    this_aa0 = aa0.loc[:,seg2_list]
    this_V = V[seg2_list]
    net_V = this_V.sum()

    this_net_aa = this_aa.copy()
    this_net_aa0 = this_aa0.copy()

    for sn in this_V.index:
        VV = this_V[sn]
        this_net_aa.loc[:,sn] = this_net_aa.loc[:,sn] * VV
        this_net_aa0.loc[:,sn] = this_net_aa0.loc[:,sn] * VV
    
    # make a Series of mean concentration in the volume
    mean_c = this_net_aa.sum(axis=1) / net_V
    mean_c0 = this_net_aa0.sum(axis=1) / net_V
    
    # find e-folding time
    td = mean_c.index.to_numpy()
    mc = mean_c.to_numpy()
    ind_ef = np.argwhere(mc < 1/np.e)[0]
    tres = td[ind_ef] # residence time in days
    
    # find e-folding time of no-reflux case
    td0 = mean_c0.index.to_numpy()
    mc0 = mean_c0.to_numpy()
    ind_ef0 = np.argwhere(mc0 < 1/np.e)[0]
    tres0 = td0[ind_ef0] # residence time in days
    
    # also calculate the "unrefluxed" residence time
    # using Qout (better because it includes rivers)
    # and doing it from the two-layer results:
    q_df = pd.read_pickle(Ldir['LOo'] + 'tef/' +
        Ldir['gtagex'] + '_' + year + '.01.01_' + year + '.12.31/' +
        'flux/two_layer_' + season + '.p')
    """
    Output: LiveOcean_output/tef/[*]/flux/two_layer_[season].p which is
        a pickled DataFrame whose index is the section names, and whose columns are:
        ['q_s', 'q_f', 'f_s', 'f_f', 's_s', 's_f', 'lon', 'lat'].
        Here _s and _f indicate that the layer is salty or fresh.
        Also we have organized all the fluxes to be positive Eastward or Northward.
    """
    if basin == 'HoodCanalInner':
        qout = -q_df.loc['hc3','q_s']
    elif basin == 'HoodCanal':
        qout = q_df.loc['hc1','q_f']
    elif basin == 'SouthSound':
        qout = q_df.loc['tn1','q_f']
    elif basin == 'PS':
        qout = -q_df.loc['ai1','q_f'] - q_df.loc['dp','q_f']
    elif basin == 'Whidbey':
        qout = -q_df.loc['wb1','q_f'] - q_df.loc['dp','q_f']
    elif basin == 'Salish':
        qout = -q_df.loc['jdf1','q_f'] + q_df.loc['sog5','q_f']
    elif basin == 'SoG':
        qout = -q_df.loc['sji1','q_f'] + q_df.loc['sog5','q_f']
    else:
        print(basin)
        print('unsupported qout')
        print('')
    tres00 = (net_V/qout)/86400
        
    
    # store results
    
    xx = x_dict[basin] + year_xdict[year]
    tres_df.loc[exp,['basin', 'year','season','tres','tres0','tres00','x']] = [basin, year, season, tres, tres0, tres00,xx]
            
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
    df.plot(x='x', y='tres0', ax=ax, marker='+',
        color=season_cdict[season], markersize=24, alpha=.7,
        linestyle='None', legend=False)

for season in season_list:
    df = tres_df[tres_df['season']==season]
    df.plot(x='x', y='tres00', ax=ax, marker='*',
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
fig.savefig(outdir + 'all_residence_times' + vftag + '.png')
plt.show()
plt.rcdefaults()

# also save the DataFrame of residence times
tres_df.pop('x') # column only needed for plotting
tres_df.to_pickle(outdir + 'tres_df' + vftag + '.p')

# and some helpful screen output
tres_basin_df = tres_df.groupby('basin').mean()
tres_basin_df['Rx Frac %'] = 100 * (tres_basin_df['tres'] - tres_basin_df['tres0']) / (tres_basin_df['tres'] - tres_basin_df['tres00'])
tres_basin_df = tres_basin_df.reindex(['Salish', 'SoG', 'PS', 'Whidbey', 'HoodCanal', 'HoodCanalInner', 'SouthSound'])
print(tres_basin_df.round())

