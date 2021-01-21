"""
Do volume, salt, and salt-squared budgets for user-specified volumes.
This is a fairly complicated task because it involves coordinating
the net salt and volume time series from flux_[get,lowpass]_s.py,
the TEF transports through any open sections of the volume,
and all the rivers flowing into the volume.

Run with a command like:
run flux_salt_budget
and then it runs for all volumes in vol_list and all years in year_list.

Currently need to do more extractions to support 2018 and 2019.

"""

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import tef_fun
import flux_fun
from importlib import reload
reload(flux_fun)

import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
from datetime import datetime, timedelta
import argparse

# debugging imports
import netCDF4 as nc
from time import time
import zfun

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', type=str, default='cas6')
parser.add_argument('-t', '--tag', type=str, default='v3')
parser.add_argument('-x', '--ex_name', type=str, default='lo8b')
parser.add_argument('-test', '--testing', type=zfun.boolean_string, default=False)
args = parser.parse_args()
testing = args.testing

old_dt = False

# Get Ldir
Ldir = Lfun.Lstart(args.gridname, args.tag)
gtagex = args.gridname + '_' + args.tag + '_' + args.ex_name

if testing:
    year_list = [2019]
    vol_list = ['Puget Sound']
    save_figs = False
else:
    year_list = [2017, 2018, 2019]
    vol_list = ['Salish Sea', 'Puget Sound', 'Hood Canal']
    save_figs = True

# make DataFrames to hold error statistics
err_df_vol = pd.DataFrame(index=year_list, columns=vol_list)
err_df_salt = pd.DataFrame(index=year_list, columns=vol_list)
err_df_salt2 = pd.DataFrame(index=year_list, columns=vol_list)

plt.close('all')
for which_vol in vol_list:

    for year in year_list:
        year_str = str(year)
    
        # select input/output location
        run_name = gtagex+'_'+year_str+'.01.01_'+year_str+'.12.31'
        indir00 = Ldir['LOo'] + 'tef2/'
        indir0 = indir00 + run_name + '/'

        outdir = indir00 + 'salt_budget_plots/'
        Lfun.make_dir(outdir)

        # load low passed segment volume, net salt, and other DataFrames
        v_lp_df = pd.read_pickle(indir0 + 'flux/daily_segment_volume.p')
        sv_lp_df = pd.read_pickle(indir0 + 'flux/daily_segment_net_salt.p')
        s2v_lp_df = pd.read_pickle(indir0 + 'flux/daily_segment_net_salt2.p')
        mix_lp_df = pd.read_pickle(indir0 + 'flux/daily_segment_mix.p')
        hmix_lp_df = pd.read_pickle(indir0 + 'flux/daily_segment_hmix.p')
        
        # new
        vt_lp_df = pd.read_pickle(indir0 + 'flux/daily_segment_vt.p')
        svt_lp_df = pd.read_pickle(indir0 + 'flux/daily_segment_svt.p')
        s2vt_lp_df = pd.read_pickle(indir0 + 'flux/daily_segment_s2vt.p')
        

        # Info specific to each volume
        # The sign for each section indicates which direction is INTO the volume.
        if which_vol == 'Salish Sea':
            seg_list = list(v_lp_df.columns)
            sect_sign_dict = {'jdf1':1, 'sog5':-1}
        elif which_vol == 'Puget Sound':
            seg_list = (flux_fun.ssA + flux_fun.ssM + flux_fun.ssT
                + flux_fun.ssS + flux_fun.ssW + flux_fun.ssH)
            sect_sign_dict = {'ai1':1, 'dp':1}
        elif which_vol == 'Hood Canal':
            seg_list = flux_fun.ssH
            sect_sign_dict = {'hc1':-1}

        v_lp_df = v_lp_df[seg_list]
        sv_lp_df = sv_lp_df[seg_list]
        s2v_lp_df = s2v_lp_df[seg_list]
        
        vt_lp_df = vt_lp_df[seg_list]
        svt_lp_df = svt_lp_df[seg_list]
        s2vt_lp_df = s2vt_lp_df[seg_list]
        
        mix_lp_df = mix_lp_df[seg_list]
        hmix_lp_df = hmix_lp_df[seg_list]
    
        sect_df = tef_fun.get_sect_df()

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
            sect_sign = sect_sign_dict[sn]
            tef_df_dict[sn], in_sign = flux_fun.get_fluxes(indir0, sn)
            if in_sign != sect_sign:
                print('WARNING: potential sign error!!')
                
        # BUDGETS
        indall = sv_lp_df.index
        
        # miscellaneous stuff that doesn't go in budget
        misc_df = pd.DataFrame(0,index=indall, columns = ['Qprism','Ftide'])
        # we initialize only columns where we sum over TEF sections
        for sect_name in tef_df_dict.keys():
            df = tef_df_dict[sect_name]
            misc_df['Qprism'] = misc_df['Qprism'] + df['Qtide']/2 # [m3/s]
            misc_df['Ftide'] = misc_df['Ftide'] + df['Ftide'] # [Watts]
        v = v_lp_df.sum(axis=1).to_numpy()
        vt = vt_lp_df.sum(axis=1).to_numpy()
        svt = svt_lp_df.sum(axis=1).to_numpy()
        s2vt = s2vt_lp_df.sum(axis=1).to_numpy()
        misc_df['V'] = v # [m3]
        
        
        # also create, from the hourly files, <s> and d/dt (int <s>2 dV)
        vh_df = pd.read_pickle(indir0 + 'flux/hourly_segment_volume.p')
        sh_df = pd.read_pickle(indir0 + 'flux/hourly_segment_salinity.p')
        vh_df = vh_df[seg_list]
        sh_df = sh_df[seg_list]
        svh_df = sh_df * vh_df # net salt in each segment (hourly)
        smeanh = svh_df.sum(axis=1).to_numpy() / vh_df.sum(axis=1).to_numpy()
        smean2h = smeanh * smeanh
        smean2vh = smean2h * vh_df.sum(axis=1).to_numpy()
        dsmean2vh_dt = np.nan + np.ones_like(smean2vh)
        dsmean2vh_dt[1:-1] = (smean2vh[2:] - smean2vh[:-2])/(2*3600)
        dsmean2vh_dt_lp = zfun.filt_godin(dsmean2vh_dt)
        pad = 36
        dsmean2v_dt_lp = dsmean2vh_dt_lp[pad:-(pad+1):24]
        
        if old_dt:
            misc_df['Smean'] = sv_lp_df.sum(axis=1)/misc_df['V']
        else:
            smeanh_lp = zfun.filt_godin(smeanh)
            smean_lp = smeanh_lp[pad:-(pad+1):24]
            misc_df['Smean'] = smean_lp
        
        # volume budget
        vol_df = pd.DataFrame(0, index=indall, columns=['Qin','Qout'])
        for sect_name in tef_df_dict.keys():
            df = tef_df_dict[sect_name]
            vol_df['Qin'] = vol_df['Qin'] + df['Qin']
            vol_df['Qout'] = vol_df['Qout'] + df['Qout']
        vol_df['Qr'] = riv_df.sum(axis=1)
        if old_dt:
            vol_df.loc[1:-1, 'dV_dt'] = (v[2:] - v[:-2]) / (2*86400)
        else:
            vol_df.loc[:, 'dV_dt'] = vt
            
        vol_df['Error'] = vol_df['dV_dt'] - vol_df.loc[:,'Qin'] - vol_df.loc[:,'Qout'] - vol_df.loc[:,'Qr']
        vol_rel_err = vol_df['Error'].mean()/vol_df['Qr'].mean()

        # salt budget
        salt_df = pd.DataFrame(0, index=indall, columns=['QSin','QSout'])
        for sect_name in tef_df_dict.keys():
            df = tef_df_dict[sect_name]
            salt_df['QSin'] = salt_df['QSin'] + df['QSin']
            salt_df['QSout'] = salt_df['QSout'] + df['QSout']
        sn = sv_lp_df[seg_list].sum(axis=1).values
        if old_dt:
            salt_df.loc[1:-1, 'dSnet_dt'] = (sn[2:] - sn[:-2]) / (2*86400)
        else:
            salt_df.loc[:, 'dSnet_dt'] = svt
            
        salt_df['Error'] = salt_df['dSnet_dt'] - salt_df['QSin'] - salt_df['QSout']
        salt_rel_err = salt_df['Error'].mean()/salt_df['QSin'].mean()
        
        misc_df['Sin'] = salt_df['QSin']/vol_df['Qin']
        misc_df['Sout'] = salt_df['QSout']/vol_df['Qout']
        Socn = misc_df['Sin'].max() # used in mixedness budget
        
        
        # re-express the salt budget using Qe notation
        salt_qe_df = pd.DataFrame(index=indall)
        salt_qe_df['dSnet_dt'] = salt_df['dSnet_dt']
        
        misc_df['Qe'] = (vol_df['Qin'] - vol_df['Qout'])/2
        misc_df['Qnet'] = -(vol_df['Qout'] + vol_df['Qin']) # like Qr, but accounting for dV/dt
        misc_df['DS'] = misc_df['Sin'] - misc_df['Sout']
        misc_df['Sbar'] = (misc_df['Sin'] + misc_df['Sout'])/2
        # NOTE: Sbar is bad notation here, need to coordinate with Variance Budget usage of Smean.
        
        salt_qe_df['QeDS'] = misc_df['Qe']*misc_df['DS']
        salt_qe_df['-QrSbar'] = -misc_df['Qnet']*misc_df['Sbar']
        salt_qe_df['Error'] = salt_qe_df['dSnet_dt'] - salt_qe_df['QeDS'] - salt_qe_df['-QrSbar']
        salt_qe_rel_err_qe = salt_qe_df['Error'].mean()/salt_qe_df['QeDS'].mean()
    
        # salinity-squared budget
        salt2_df = pd.DataFrame(0, index=indall, columns=['QS2in','QS2out'])
        for sect_name in tef_df_dict.keys():
            df = tef_df_dict[sect_name]
            salt2_df['QS2in'] = salt2_df['QS2in'] + df['QS2in']
            salt2_df['QS2out'] = salt2_df['QS2out'] + df['QS2out']
        s2n = s2v_lp_df[seg_list].sum(axis=1).to_numpy()
        if old_dt:
            ds2net_dt = (s2n[2:] - s2n[:-2]) / (2*86400)
            salt2_df.loc[1:-1, 'dS2net_dt'] = ds2net_dt
        else:
            ds2net_dt = s2vt
            salt2_df.loc[:, 'dS2net_dt'] = ds2net_dt
        
        salt2_df['VMix'] = mix_lp_df[seg_list].sum(axis=1).to_numpy()
        salt2_df['HMix'] = hmix_lp_df[seg_list].sum(axis=1).to_numpy()
        salt2_df['Mix_resolved'] = salt2_df['VMix'] + salt2_df['HMix']
        
        misc_df['S2in'] = salt2_df['QS2in']/vol_df['Qin']
        misc_df['S2out'] = salt2_df['QS2out']/vol_df['Qout']
        misc_df['S2mean'] = s2v_lp_df.sum(axis=1)/misc_df['V']
    
        salt2_df['Mix_numerical'] = salt2_df['dS2net_dt'] - salt2_df['QS2in'] - salt2_df['QS2out'] - salt2_df['Mix_resolved']
        salt2_df['Mixing'] = salt2_df['Mix_numerical'] + salt2_df['Mix_resolved']
        salt2_rel_err = salt2_df['Mix_numerical'].mean()/salt2_df['Mix_resolved'].mean()
        
        # variance budget for S'^2
        sp2_df = pd.DataFrame(index=indall)
        
        misc_df['Sp2in'] = misc_df['S2in'] - 2*misc_df['Smean']*misc_df['Sin'] + misc_df['Smean']**2
        misc_df['Sp2out'] = misc_df['S2out'] - 2*misc_df['Smean']*misc_df['Sout'] + misc_df['Smean']**2
        
        sp2_df['QSp2in'] = misc_df['Sp2in'] * vol_df['Qin']
        sp2_df['QSp2out'] = misc_df['Sp2out'] * vol_df['Qout']
        sp2_df['QrSp2'] = vol_df['Qr']*misc_df['S2mean']
    
        if old_dt:
            smean2_net = misc_df['Smean'].to_numpy()**2 * misc_df['V'].to_numpy()
            dsmean2_net_dt = (smean2_net[2:] - smean2_net[:-2]) / (2*86400)
        else:
            dsmean2_net_dt = dsmean2v_dt_lp
            
        if old_dt:
            sp2_df.loc[1:-1,'dSp2net_dt'] = ds2net_dt - dsmean2_net_dt
        else:
            sp2_df.loc[:,'dSp2net_dt'] = ds2net_dt - dsmean2_net_dt
            
        sp2_df['Mix_resolved'] = salt2_df['Mix_resolved']
        sp2_df['Mix_numerical'] = sp2_df['dSp2net_dt'] - sp2_df['QSp2in'] - sp2_df['QSp2out'] - sp2_df['QrSp2'] - sp2_df['Mix_resolved']
        sp2_df['Mixing'] = sp2_df['Mix_resolved'] + sp2_df['Mix_numerical']
        
        # Normalized S'^2 budget
        QSS = vol_df['Qr'].mean() * misc_df['Sin'] * misc_df['Sin']
        sp2n_df = sp2_df.copy()
        for cn in sp2n_df.columns:
            sp2n_df[cn] = sp2n_df[cn]/QSS
            
        # Mixedness budget
        m_df = pd.DataFrame(index=indall)
        m_df['dmnet_dt'] = Socn*salt_df['dSnet_dt'] - salt2_df['dS2net_dt']
        m_df['Qmin'] = vol_df['Qin']*(Socn*misc_df['Sin'] - misc_df['S2in'])
        m_df['Qmout'] = vol_df['Qout']*(Socn*misc_df['Sout'] - misc_df['S2out'])
        m_df['Mixing'] = -salt2_df['Mixing']
        m_df['Error'] = m_df['dmnet_dt'] - m_df['Qmin'] - m_df['Qmout'] - m_df['Mixing']
        
        # normalized mixedness
        mn_df = m_df * (365 * 86400) / (misc_df['V'].mean() * (Socn/2) * (Socn/2))
        
        # make sure everything is numeric
        for cn in misc_df.columns:
            misc_df[cn] = pd.to_numeric(misc_df[cn])
        for cn in vol_df.columns:
            vol_df[cn] = pd.to_numeric(vol_df[cn])
        for cn in salt_df.columns:
            salt_df[cn] = pd.to_numeric(salt_df[cn])
        for cn in salt2_df.columns:
            salt2_df[cn] = pd.to_numeric(salt2_df[cn])
        for cn in salt_qe_df.columns:
            salt_qe_df[cn] = pd.to_numeric(salt_qe_df[cn])
        for cn in sp2_df.columns:
            sp2_df[cn] = pd.to_numeric(sp2_df[cn])
        for cn in sp2n_df.columns:
            sp2n_df[cn] = pd.to_numeric(sp2n_df[cn])
        for cn in m_df.columns:
            m_df[cn] = pd.to_numeric(m_df[cn])
        for cn in m_df.columns:
            mn_df[cn] = pd.to_numeric(mn_df[cn])
            
        # save DataFrames to pickle files
        misc_df.to_pickle(outdir + 'misc_df_' + year_str + '_' + which_vol.replace(' ','_') + '.p')
        vol_df.to_pickle(outdir + 'vol_df_' + year_str + '_' + which_vol.replace(' ','_') + '.p')
        salt_df.to_pickle(outdir + 'salt_df_' + year_str + '_' + which_vol.replace(' ','_') + '.p')
        salt2_df.to_pickle(outdir + 'salt2_df_' + year_str + '_' + which_vol.replace(' ','_') + '.p')
        sp2_df.to_pickle(outdir + 'sp2_df_' + year_str + '_' + which_vol.replace(' ','_') + '.p')
        sp2n_df.to_pickle(outdir + 'sp2n_df_' + year_str + '_' + which_vol.replace(' ','_') + '.p')
        m_df.to_pickle(outdir + 'm_df_' + year_str + '_' + which_vol.replace(' ','_') + '.p')
            
        # accumulate error statistics (all %)
        err_df_vol.loc[year, which_vol] = vol_rel_err*100
        err_df_salt.loc[year, which_vol] = salt_rel_err*100
        err_df_salt2.loc[year, which_vol] = salt2_rel_err*100
            
        # debugging: get net salt and salt-squared flux as a check
        if False:
            salt_df['QSnet'] = 0
            salt2_df['QS2net'] = 0
            tt0 = time()
            for sect_name in sect_sign_dict.keys():
                print('Debugging - processing ' + sect_name)
                sys.stdout.flush()
                sect_sign = sect_sign_dict[sect_name]
                ext_fn = indir0 + 'extractions/' + sect_name + '.nc'
                ext = nc.Dataset(ext_fn)
                ext_q = sect_sign * ext['q'][:]
                ext_salt = ext['salt'][:]
                xqs = ((ext_q * ext_salt).sum(axis=2)).sum(axis=1)
                xqs2 = ((ext_q * ext_salt * ext_salt).sum(axis=2)).sum(axis=1)
                xqs_lp = zfun.filt_godin(xqs)
                xqs2_lp = zfun.filt_godin(xqs2)
                pad = 36
                xqs_lp = xqs_lp[pad:-(pad+1):24]
                xqs2_lp = xqs2_lp[pad:-(pad+1):24]
                salt_df['QSnet'] += xqs_lp
                salt2_df['QS2net'] += xqs2_lp
                print('  -- took %0.2f sec' % (time()-tt0))
                sys.stdout.flush()
            salt_df['Error_alt'] = salt_df['dSnet_dt'] - salt_df['QSnet']
            salt2_df['Mix_numerical_alt'] = salt2_df['dS2net_dt'] - salt2_df['QS2net'] - salt2_df['Mix_resolved']
            
        # --------------- Volume and Salt Budgets --------------------------------------
        fig1 = plt.figure(figsize=(16,7))
        
        ax = fig1.add_subplot(121)
        vol_df[['dV_dt','Qin','Qout', 'Qr','Error']].plot(ax=ax, grid=True).legend(loc='upper right')
        ax.set_title('Volume Budget (m3/s)')
        ax.text(.05,.9, 'Mean Error / Mean Qr = %0.2f%%' % (vol_rel_err*100), transform=ax.transAxes, fontsize=14)

        ax = fig1.add_subplot(122)
        #salt_df[['dSnet_dt','QSin','QSout','Error']].plot(ax=ax, grid=True).legend(loc='upper right')
        salt_df.plot(ax=ax, grid=True).legend(loc='upper right')
        ax.set_title(year_str + ' ' + which_vol + ' Salt Budget (g/kg m3/s)')
        ax.text(.05,.9, 'Mean Error / Mean QSin = %0.2f%%' % (salt_rel_err*100), transform=ax.transAxes, fontsize=14)
        
        if save_figs:
            outname1 = outdir + 'vol_salt_budget_' + year_str + '_' + which_vol.replace(' ','_')
            fig1.savefig(outname1 + '.png')
        
        # --------------- Salinity-squared Budget ------------------------------------------
        fig2 = plt.figure(figsize=(16,7))
        
        ax = fig2.add_subplot(121)
        #salt2_df[['dS2net_dt','QS2in','QS2out','Mix_resolved','Mix_numerical']].plot(ax=ax, grid=True).legend(loc='upper right')
        salt2_df.plot(ax=ax, grid=True).legend(loc='upper right')
        ax.set_title(year_str + ' ' + which_vol + ' Salt2 Budget (g2/kg2 m3/s)')
                
        ax = fig2.add_subplot(122)
        salt2_df[['Mix_resolved','VMix','HMix','Mix_numerical']].plot(ax=ax, grid=True).legend(loc='upper right')
        ax.set_title(year_str + ' ' + which_vol + ' Salt2 Budget (g2/kg2 m3/s)')
        ax.text(.05,.9, 'Mean Error / Mean Mix_resolved = %0.2f%%' % (salt2_rel_err*100), transform=ax.transAxes, fontsize=14)

        if save_figs:
            outname2 = outdir + 'salt2_budget_' + year_str + '_' + which_vol.replace(' ','_')
            fig2.savefig(outname2 + '.png')

        # --------------- S'^2 squared Budget ----------------------------------------
        fig3 = plt.figure(figsize=(12,7))
                
        ax = fig3.add_subplot(111)
        sp2n_df[['dSp2net_dt', 'QSp2in', 'QSp2out', 'QrSp2', 'Mixing', 'Mix_resolved']].plot(ax=ax, grid=True).legend(loc='upper right')
        ax.set_title(year_str + ' ' + which_vol + ' Sp2 Budget Normalized')
        ax.set_ylim(-4,4)
        
        if save_figs:
            outname3 = outdir + 'sp2n_budget_' + year_str + '_' + which_vol.replace(' ','_')
            fig3.savefig(outname3 + '.png')
        
        
        # --------------- Mixedness Budget ------------------------------------
        fig4 = plt.figure(figsize=(12,7))
                
        ax = fig4.add_subplot(211)
        m_df[['dmnet_dt', 'Qmin', 'Qmout', 'Mixing', 'Error']].plot(ax=ax, grid=True).legend(loc='upper right')
        ax.set_title(year_str + ' ' + which_vol + ' Mixedness Budget')
        
        ax = fig4.add_subplot(212)
        misc_df['Qprism'].plot(ax=ax, grid=True)
        #ax.set_ylim(-4,4)
        
        if save_figs:
            outname4 = outdir + 'm_budget_' + year_str + '_' + which_vol.replace(' ','_')
            fig4.savefig(outname4 + '.png')

        # --------------- Normalized Mixedness Budget ------------------------------------
        fig5 = plt.figure(figsize=(12,7))
                
        ax = fig5.add_subplot(211)
        mn_df[['dmnet_dt', 'Qmin', 'Qmout', 'Mixing', 'Error']].plot(ax=ax, grid=True).legend(loc='upper right')
        ax.set_title(year_str + ' ' + which_vol + ' Normalized Mixedness Budget')
        ax.set_ylabel('V*(Socn/2)^2 per year')
        
        ax = fig5.add_subplot(212)
        misc_df['Qprism'].plot(ax=ax, grid=True)
        #ax.set_ylim(-4,4)
        
        if save_figs:
            outname5 = outdir + 'm_budget_' + year_str + '_' + which_vol.replace(' ','_')
            fig5.savefig(outname5 + '.png')
        
        
        # --------------- Checking ---------------------------------------------------
        fig6 = plt.figure(figsize=(12,6))

        ax = fig6.add_subplot(111)
        salt2_df[['Mixing']].plot(ax=ax, grid=True, label='from salt2',style='-r').legend(loc='upper right')
        sp2_df[['Mixing']].plot(ax=ax, grid=True, label='from sp2',style='-b').legend(loc='upper right')
        ax.set_title(year_str + ' ' + which_vol + ' Mixing (g2/kg2 m3/s)')
        ax.text(.05, .9, 'Mean Mixing from salt2 = %.5e' % (salt2_df['Mixing'].mean()), transform=ax.transAxes, c='r')
        ax.text(.05, .8, 'Mean Mixing from sp2 = %.5e' % (sp2_df['Mixing'].mean()), transform=ax.transAxes, c='b')
        ax.text(.05, .2, 'Ratio = %.3f' % (salt2_df['Mixing'].mean()/sp2_df['Mixing'].mean()), transform=ax.transAxes, c='purple')
        
plt.show()

