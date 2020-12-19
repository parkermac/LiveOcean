"""
Make low-passed versions of the hourly extractions made by
flux_get_s.py.

Takes about 30 sec per year.

"""

# imports
import pandas as pd
import numpy as np

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
Ldir = Lfun.Lstart('cas6', 'v3')

indir0 = Ldir['LOo'] + 'tef2/'

for year in [2017, 2018, 2019]:
    print('Working on ' + str(year))
    
    outname = 'cas6_v3_lo8b_' + str(year) + '.01.01_' + str(year) + '.12.31'
    outdir = indir0 + outname + '/flux/'

    # load DataFrames of hourly volume and salinity
    v_df = pd.read_pickle(outdir + 'hourly_segment_volume.p')
    s_df = pd.read_pickle(outdir + 'hourly_segment_salinity.p')
    mix_df = pd.read_pickle(outdir + 'hourly_segment_mix.p')
    hmix_df = pd.read_pickle(outdir + 'hourly_segment_hmix.p')
    s2_df = pd.read_pickle(outdir + 'hourly_segment_salinity2.p')
    seg_list = list(v_df.columns)
    
    # Form tidally averaged time series of net salt in segments.
    # We do this exactly the same as
    # it was done in bulk_calc.py so the resulting time indices are identical.
    pad = 36
    dt_list = list(v_df.index)
    dt_list = dt_list[pad:-(pad+1):24]
    s_lp_df = pd.DataFrame(index=dt_list, columns=seg_list)
    v_lp_df = pd.DataFrame(index=dt_list, columns=seg_list)
    sv_lp_df = pd.DataFrame(index=dt_list, columns=seg_list)
    s2_lp_df = pd.DataFrame(index=dt_list, columns=seg_list)

    vt_lp_df = pd.DataFrame(index=dt_list, columns=seg_list)
    svt_lp_df = pd.DataFrame(index=dt_list, columns=seg_list)
    s2vt_lp_df = pd.DataFrame(index=dt_list, columns=seg_list)

    mix_lp_df = pd.DataFrame(index=dt_list, columns=seg_list)
    hmix_lp_df = pd.DataFrame(index=dt_list, columns=seg_list)
    s2v_lp_df = pd.DataFrame(index=dt_list, columns=seg_list)
    for seg_name in seg_list:
        
        v = v_df.loc[:,seg_name].values
        s = s_df.loc[:,seg_name].values
        sv = s*v
        
        mix = mix_df.loc[:,seg_name].values
        hmix = hmix_df.loc[:,seg_name].values
        s2 = s2_df.loc[:,seg_name].values
        s2v = s2*v
        
        # new 2020.12.18 form time derivatives
        vt = np.nan + np.ones_like(v)
        vt[1:-1] = (v[2:] - v[:-2])/(2*3600)
        svt = np.nan + np.ones_like(sv)
        svt[1:-1] = (sv[2:] - sv[:-2])/(2*3600)
        s2vt = np.nan + np.ones_like(s2v)
        s2vt[1:-1] = (s2v[2:] - s2v[:-2])/(2*3600)
        
        v_lp = zfun.filt_godin(v)
        s_lp = zfun.filt_godin(s)
        sv_lp = zfun.filt_godin(sv)
        
        vt_lp = zfun.filt_godin(vt)
        svt_lp = zfun.filt_godin(svt)
        s2vt_lp = zfun.filt_godin(s2vt)
        
        mix_lp = zfun.filt_godin(mix)
        hmix_lp = zfun.filt_godin(hmix)
        s2_lp = zfun.filt_godin(s2)
        s2v_lp = zfun.filt_godin(s2v)
        
        # subsample and cut off nans
        v_lp = v_lp[pad:-(pad+1):24]
        s_lp = s_lp[pad:-(pad+1):24]
        sv_lp = sv_lp[pad:-(pad+1):24]
        
        vt_lp = vt_lp[pad:-(pad+1):24]
        svt_lp = svt_lp[pad:-(pad+1):24]
        s2vt_lp = s2vt_lp[pad:-(pad+1):24]
        
        mix_lp = mix_lp[pad:-(pad+1):24]
        hmix_lp = hmix_lp[pad:-(pad+1):24]
        s2_lp = s2_lp[pad:-(pad+1):24]
        s2v_lp = s2v_lp[pad:-(pad+1):24]
        
        # save to DataFrames
        v_lp_df.loc[:,seg_name] = v_lp
        s_lp_df.loc[:,seg_name] = s_lp
        sv_lp_df.loc[:,seg_name] = sv_lp
        
        vt_lp_df.loc[:,seg_name] = vt_lp
        svt_lp_df.loc[:,seg_name] = svt_lp
        s2vt_lp_df.loc[:,seg_name] = s2vt_lp
        
        mix_lp_df.loc[:,seg_name] = mix_lp
        hmix_lp_df.loc[:,seg_name] = hmix_lp
        s2_lp_df.loc[:,seg_name] = s2_lp
        s2v_lp_df.loc[:,seg_name] = s2v_lp
    
    # save results to disk.
    v_lp_df.to_pickle(outdir + 'daily_segment_volume.p')
    s_lp_df.to_pickle(outdir + 'daily_segment_salinity.p')
    sv_lp_df.to_pickle(outdir + 'daily_segment_net_salt.p')
    
    vt_lp_df.to_pickle(outdir + 'daily_segment_vt.p')
    svt_lp_df.to_pickle(outdir + 'daily_segment_svt.p')
    s2vt_lp_df.to_pickle(outdir + 'daily_segment_s2vt.p')
    
    mix_lp_df.to_pickle(outdir + 'daily_segment_mix.p')
    hmix_lp_df.to_pickle(outdir + 'daily_segment_hmix.p')
    s2_lp_df.to_pickle(outdir + 'daily_segment_salinity2.p')
    s2v_lp_df.to_pickle(outdir + 'daily_segment_net_salt2.p')
    
