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
#import netCDF4 as nc
import argparse

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
import zrfun

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

import tef_fun
import flux_fun

from importlib import reload
reload(flux_fun)

def boolean_string(s):
    if s not in ['False', 'True']:
        raise ValueError('Not a valid boolean string')
    return s == 'True' # note use of ==

# optional command line arguments, can be input in any order
parser = argparse.ArgumentParser()
parser.add_argument('-src', '--source', nargs='?', type=str, default='ocean')
parser.add_argument('-sink', '--sinking', default=False, type=boolean_string)
parser.add_argument('-conv', '--check_convergence', default=False, type=boolean_string)
args = parser.parse_args()

source = args.source
sinking = args.sinking
check_convergence = args.check_convergence

print('Running integration with source = ' + source)
# Valid choices for source:
#
# Dye from the ocean meant to reproduce the actual mean salinity
# - ocean_salt
#
# Dye coming in with concentration = 1:
# - ocean
# - river_fraser
# - river_deschutes
# - river_skagit
#
# Initial condition, dye concentration = 1 in a basin
# - ic_hood_canal_inner
# - MORE: see flux_fun.ic_seg2_dict
#

# Input directory
indir0 = Ldir['LOo'] + 'tef/'
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/flux/'

# hacky way of getting the year, assumes "item" is of the form:
# 'cas6_v3_lo8b_2017.01.01_2017.12.31'
year_str = item.split('_')[-1].split('.')[0]
year = int(year_str)

# set ocean salinity BC by year
if year == 2017:
    os_jdf = 33.1; os_sog = 30.5
elif year == 2018:
    os_jdf = 33.18; os_sog = 30.8
elif year == 2019:
    os_jdf = 33.34; os_sog = 31.0

plt.close('all')

# loop over all seasons
for season in flux_fun.season_list:
    
    print('===== ' + season + ' =====')

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
        if source == 'ocean_salt':
            if 'J1' in seg_name:
                f.loc[seg_name,'ocean_s'] = os_jdf
            elif 'G6' in seg_name:
                f.loc[seg_name,'ocean_s'] = os_sog
        elif source == 'ocean':
            if 'J1' in seg_name:
                f.loc[seg_name,'ocean_s'] = 1
            elif 'G6' in seg_name:
                f.loc[seg_name,'ocean_s'] = os_sog/os_jdf
        elif source == 'river_fraser':
            if 'G3' in seg_name:
                f.loc[seg_name,'river_f'] = 1
        elif source == 'river_deschutes':
            if 'S4' in seg_name:
                f.loc[seg_name,'river_f'] = 1
        elif source == 'river_skagit':
            if 'W4' in seg_name:
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

    if 'ocean' in source:
        ff[:,0] = f.loc[:,'ocean_s'].values # column 0 is the ocean inflow
    elif 'river' in source:
        ff[:,3] = f.loc[:,'river_f'].values # column 3 is the river inflow
    
    if 'ic' in source:
        seg2_list = flux_fun.ic_seg2_dict[source]
        for seg_name in f.index:
            if seg_name in seg2_list:
            # if source in ['ic_hood_canal', 'ic_hood_canal']:
            #     this_seg_list = ['H'+ str(item) for item in range(3,9)]
            #     if seg_name[:2] in this_seg_list:
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
        # here we are updating the concentration array in the segment bins
        # but not in the ocean and river bins
        
        # NOTE: do this with a flag on the command line instead of ic0 name?
        # here is where will intervene to do a different sort of "ic" experiment
        # to zero-out all concentrations outside of the release bins
        # bin_list = []
        # if 'ic0' in source:
        #     for seg_name in this_seg_list:
        #         bin_list.append(seg_name + '_s')
        #         bin_list.append(seg_name + '_f')
        #     bin_list_full = list(q_df.index)
        #     ic_filter = np.zeros(NR)
        #     # this is a filter of zeros except it it ones where we set the IC
        #     for bin in bin_list_full:
        #         if bin in bin_list:
        #             ibin = bin_list_full.index(bin)
        #             ic_filter[ibin] = 1
        #     # and finally we zero out all the concentrations outside of the IC region
        #     ff[:,4:] = np.tile(ic_filter*c,(NR,1))
        #     # this overrides the update above
            
        
        # similar calculations to the above but with a tracer that increases
        # over time at a rate that appears to be 100% per year
        qffa = (q*ffa).sum(axis=1)
        ca = ca + dt*ivv*qffa + dt*c/(365*86400)
        ffa[:,4:] = np.tile(ca,(NR,1))
        
    if check_convergence:
        fig = plt.figure(figsize=(12,8))
        ax = fig.add_subplot(111)
        ax.plot(t_arr/(365*86400), c_arr, '-')
        ax.set_xlabel('Time (years)')
        ax.set_xlim(0,nyears)
        ax.set_title(season.title())
        plt.show()

    # create and fill the output DataFrame
    #
    # this one is just the final state (useful for age)
    cc = pd.DataFrame(index=q_df.index,columns=['c', 'ca'])
    cc['c'] = c
    cc['ca'] = ca

    # and another one for the time-dependent fields in an array aa
    aa = pd.DataFrame(c_arr, index=t_arr/86400, columns=q_df.index)
    # note that the index is time (days) and the columns are the bin names

    sink_tag = ''
    if sinking:
        sink_tag = '_sinking'
        
    cc.to_pickle(indir + 'cc_' + source + '_' + season + sink_tag + '.p')
    aa.to_pickle(indir + 'aa_' + source + '_' + season + sink_tag + '.p')
    
    