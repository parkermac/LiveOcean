"""
Run the flux engine!

In the comments the word "bin" is used to mean an individual volume in the
system, so there are two bins per segment, an upper one and a lower one.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
import argparse

def boolean_string(s):
    if s not in ['False', 'True']:
        raise ValueError('Not a valid boolean string')
    return s == 'True' # note use of ==

# optional command line arguments, can be input in any order
parser = argparse.ArgumentParser()

# this one is essential
parser.add_argument('-src', '--source', nargs='?', type=str, default='')

# these are optional
parser.add_argument('-sink', '--sinking', default=False, type=boolean_string)

# typically don't need to change these because I am analyzing one LO instance
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas6')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v3')
parser.add_argument('-x', '--exname', nargs='?', type=str, default='lo8b')
args = parser.parse_args()

if len(args.source) == 0:
    print('** Need to specify a source at the command line **')
    sys.exit()

source = args.source
sinking = args.sinking

# get Ldir and associated LO modules
import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname=args.gridname, tag=args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.exname
import zrfun

import tef_fun
import flux_fun
from importlib import reload
reload(flux_fun)

print('Running integration with source = ' + source)
# Valid choices for source:
#
# Dye from the ocean meant to reproduce the actual mean salinity
# - S_OceanSalt
#
# Dye coming in with concentration = 1:
# - S_Ocean
# - S_FraserRiver
# - S_DeschutesRiver
# - S_SkagitRiver
# - S_AllRiver
#
# Initial condition, dye concentration = 1 in part of a basin
# - IC_HoodCanalInner
#
# - MORE: see flux_fun.ic_seg2_dict

if source == 'S_OceanSalt':
    # this one is meant only for validation, and the validation is only
    # meaningful when we are averaging over a time when the mean state of
    # the system is approximately cyclic, not steadily increasing or
    # decreasing like it may for a singel season.
    season_list = ['full']
else:
    season_list = ['winter','spring','summer','fall']

# Input directory
indir0 = Ldir['LOo'] + 'tef/'

# Output directory
outdir0 = indir0 + 'flux_engine/'
Lfun.make_dir(outdir0)
outdir = outdir0 + Ldir['gtagex'] + '/'
Lfun.make_dir(outdir)

plt.close('all')
for year in [2017, 2018, 2019]:
    
    date_range = str(year) + '.01.01_' + str(year) + '.12.31'
    extraction_name = Ldir['gtagex'] + '_' + date_range
    indir = indir0 + extraction_name + '/flux/'
    
    # set ocean salinity BC by year
    if year == 2017:
        os_jdf = 33.1; os_sog = 30.5
    elif year == 2018:
        os_jdf = 33.18; os_sog = 30.8
    elif year == 2019:
        os_jdf = 33.34; os_sog = 31.0

    # loop over all seasons
    for season in season_list:
            
        # form the core of the output name
        if sinking == True:
            source_str = source + 'Sink'
        else:
            source_str = source
        outname_main = source_str + '_' + str(year) + '_' + season
        print('=== ' + outname_main + ' ===')
        
        # load DateFrames of transport and volume
        q_df = pd.read_pickle(indir + 'q_df_' + season + '.p')
        v_df = pd.read_pickle(indir + 'volumes.p')
        V = flux_fun.get_V(v_df)
    
        # "f" is a DataFrame organized like q_df but whose entries
        # are the forced values of the tracer
        f = pd.DataFrame(0, index=q_df.index, columns=q_df.columns)
        # initialized here with zeros everywhere?

        # set forcing values of ocean or river boundary conditions
        for seg_name in f.index:
            if source == 'S_OceanSalt':
                if 'J1' in seg_name:
                    f.loc[seg_name,'ocean_s'] = os_jdf
                elif 'G6' in seg_name:
                    f.loc[seg_name,'ocean_s'] = os_sog
            elif source == 'S_Ocean':
                if 'J1' in seg_name:
                    f.loc[seg_name,'ocean_s'] = 1
                elif 'G6' in seg_name:
                    f.loc[seg_name,'ocean_s'] = 1 #os_sog/os_jdf
            elif source == 'S_FraserRiver':
                if 'G3' in seg_name:
                    f.loc[seg_name,'river_f'] = 1
            elif source == 'S_DeschutesRiver':
                if 'S4' in seg_name:
                    f.loc[seg_name,'river_f'] = 1
            elif source == 'S_SkagitRiver':
                if 'W4_f' in seg_name:
                    f.loc[seg_name,'river_f'] = 1
            elif source == 'S_AllRiver':
                if '_f' in seg_name:
                    f.loc[seg_name,'river_f'] = 1

        # we will do the integration just with numpy arrays
        NR = len(q_df.index) # equal to the number of bins (2 x number of segments)
        NC = len(q_df.columns) # same as NR but plus 4 for ocean and river open boundaries
        vv = V.values
        ivv = 1/vv

        c = np.zeros(NR) # the final concentration in all bins
        ca = np.zeros(NR) # the final concentration in all bins for an "aging" tracer
        q = q_df.values # transports (m3/s) ROW=to, COLUMNS=from
        ff = np.zeros((NR,NC)) # how is ff different from f (ff=array, f=DataFrame)?
        ffa = np.zeros((NR,NC))

        if 'Ocean' in source:
            ff[:,0] = f.loc[:,'ocean_s'].values # column 0 is the ocean inflow
        elif 'River' in source:
            ff[:,3] = f.loc[:,'river_f'].values # column 3 is the river inflow
    
        if 'IC_' in source:
            seg2_list = flux_fun.ic_seg2_dict[source]
            for seg_name in f.index:
                if seg_name in seg2_list:
                    jj = int(np.argwhere(f.index==seg_name))
                    ff[jj,jj+4] = 1
                    c[jj] = 1

        dt = 3600 #3e3 # time step (seconds)
        nyears = 6
        NT = int(nyears*365*86400/dt) # number of time steps
        savedays = 2#10 # number of days between saves
        Nsave = int(savedays*86400/dt)
        # Nsave = number of time steps between saves in order to save every "savedays"

        # arrays to STORE time-varying information
        c_arr = np.nan + np.ones((int(NT/Nsave)+1, NR))
        t_arr = np.nan + np.ones(int(NT/Nsave)+1)
    
        dz = .5e-4 #1e-4 # a parameter to control sinking rate

        for ii in range(NT):
    
            # save the bin concentration to an array (and save time) every "savedays" days
            if np.mod(ii,Nsave)==0 :
                c_arr[int(ii/Nsave), :] = c.copy()
                t_arr[int(ii/Nsave)] = dt * ii
            
            # General scheme: q is the fluxes and ff is the concentration at this
            # time step.
            qff = (q*ff).sum(axis=1)
            # So this sum of transport (q) times concentration (ff) across all the columns
            # gives the transport rate into a given bin
        
            c = c + dt*ivv*qff
            # Here we advance the concentration using dc/dt = net inflow / volume
        
            if sinking == True:
                NC2 = int(len(c)/2)
                for jj in range(NC2):
                    c_f = c[2*jj + 1]
                    c[2*jj + 1] -= c_f*dz
                    c[2*jj] += c_f*dz
            ff[:,4:] = np.tile(c,(NR,1))
                    
            # similar calculations to the above but with a tracer that increases
            # over time at a rate that appears to be 100% per year
            qffa = (q*ffa).sum(axis=1)
            ca = ca + dt*ivv*qffa + dt*c/(365*86400)
            ffa[:,4:] = np.tile(ca,(NR,1))
        

        # create and fill the output DataFrames
        
        if 'IC_' not in source:
            # this one is just the final state (useful for age)
            cc = pd.DataFrame(index=q_df.index,columns=['c', 'ca'])
            cc['c'] = c
            cc['ca'] = ca
            cc.to_pickle(outdir + outname_main + '_AGE.p')

        # and another one for the time-dependent fields in an array aa
        aa = pd.DataFrame(c_arr, index=t_arr/86400, columns=q_df.index)
        # note that the index is time (days) and the columns are the bin names
        aa.to_pickle(outdir + outname_main + '.p')
    
    