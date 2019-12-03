"""
Make low-passed version of the hourly extraction made by
flux_get_s.py.

"""

# imports
import pandas as pd

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
Ldir = Lfun.Lstart('cas6', 'v3')


# select input/output location
indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

# load DataFrames of hourly volume and salinity
v_df = pd.read_pickle(indir + 'hourly_segment_volume.p')
s_df = pd.read_pickle(indir + 'hourly_segment_salinity.p')
seg_list = list(v_df.columns)
    
# Form tidally averaged time series of net salt in segments.  We do this exactly the same as
# it was done in bulk_calc.py so the resulting time indices are identical.
pad = 36
dt_list = list(v_df.index)
dt_list = dt_list[pad:-(pad+1):24]
s_lp_df = pd.DataFrame(index=dt_list, columns=seg_list)
v_lp_df = pd.DataFrame(index=dt_list, columns=seg_list)
sv_lp_df = pd.DataFrame(index=dt_list, columns=seg_list)
for seg_name in seg_list:
    v = v_df.loc[:,seg_name].values
    s = s_df.loc[:,seg_name].values
    sv = s*v
    v_lp = zfun.filt_godin(v)
    s_lp = zfun.filt_godin(s)
    sv_lp = zfun.filt_godin(sv)
    # subsample and cut off nans
    v_lp = v_lp[pad:-(pad+1):24]
    s_lp = s_lp[pad:-(pad+1):24]
    sv_lp = sv_lp[pad:-(pad+1):24]
    # save to DataFrames
    s_lp_df.loc[:,seg_name] = s_lp
    v_lp_df.loc[:,seg_name] = v_lp
    sv_lp_df.loc[:,seg_name] = sv_lp
    
# save results to disk.
s_lp_df.to_pickle(indir + 'daily_segment_salinity.p')
v_lp_df.to_pickle(indir + 'daily_segment_volume.p')
sv_lp_df.to_pickle(indir + 'daily_segment_net_salt.p')
