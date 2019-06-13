# -*- coding: utf-8 -*-
"""
Code to calculate a TEF time series using Marvin Lorenz' new multi-layer code.

Based on his code, and modified by PM.
"""

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun
Ldir = Lfun.Lstart()

import numpy as np
import pickle
import matplotlib.pyplot as plt

import tef_fun_lorenz as tfl
from importlib import reload
reload(tfl)

# choose input and organize output
Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'tef/'
# choose the tef extraction to process
item = Lfun.choose_item(indir0)
indir0 = indir0 + item + '/'
indir = indir0 + 'processed/'
sect_list_raw = os.listdir(indir)
sect_list_raw.sort()
sect_list = [item for item in sect_list_raw if ('.p' in item)]
print(20*'=' + ' Processed Sections ' + 20*'=')
print(*sect_list, sep=", ")
print(61*'=')
# select which sections to process
my_choice = input('-- Input section to process (e.g. sog5, or Return to process all): ')
if len(my_choice)==0:
    # full list
    pass
else: # single item
    if (my_choice + '.p') in sect_list:
        sect_list = [my_choice + '.p']
    else:
        print('That section is not available')
        sys.exit()
outdir = indir0 + 'bulk/'
Lfun.make_dir(outdir)

testing = False
for snp in sect_list:
    print('Working on ' + snp)
    out_fn = outdir + snp

    # load the data file
    tef_ex=pickle.load(open(indir + snp, 'rb'))
    # Notes on the data:
    # data.keys() => dict_keys(['tef_q', 'tef_qs', 'sbins', 'ot', 'qnet', 'fnet', 'ssh'])
    # data['tef_q'].shape => (8761, 1000), so packed [hour, salinity bin]
    # sbins are packed low to high
    # ot is time in seconds from 1/1/1970
    sbins = tef_ex['sbins']
    ot = tef_ex['ot']
    tef_q = tef_ex['tef_q']
    tef_qs = tef_ex['tef_qs']
    qnet = tef_ex['qnet']
    fnet = tef_ex['fnet']
    ssh = tef_ex['ssh']

    # low-pass
    if True:
        # tidal averaging
        tef_q_lp = zfun.filt_godin_mat(tef_q)
        tef_qs_lp = zfun.filt_godin_mat(tef_qs)
        qnet_lp = zfun.filt_godin(qnet)
        fnet_lp = zfun.filt_godin(fnet)
        ssh_lp = zfun.filt_godin(ssh)
        pad = 36
    else:
        # nday Hanning window
        nday = 120
        nfilt = nday*24
        tef_q_lp = zfun.filt_hanning_mat(tef_q, n=nfilt)
        tef_qs_lp = zfun.filt_hanning_mat(tef_qs, n=nfilt)
        qnet_lp = zfun.filt_hanning(qnet, n=nfilt)
        fnet_lp = zfun.filt_hanning(fnet, n=nfilt)
        ssh_lp = zfun.filt_hanning(ssh, n=nfilt)
        pad = int(np.ceil(nfilt/2))

    # subsample and cut off nans
    tef_q_lp = tef_q_lp[pad:-(pad+1):24, :]
    tef_qs_lp = tef_qs_lp[pad:-(pad+1):24, :]
    ot = ot[pad:-(pad+1):24]
    qnet_lp = qnet_lp[pad:-(pad+1):24]
    fnet_lp = fnet_lp[pad:-(pad+1):24]
    ssh_lp = ssh_lp[pad:-(pad+1):24]

    # get sizes and make sedges (the edges of sbins)
    DS=sbins[1]-sbins[0]
    sedges = np.concatenate((sbins,np.array([sbins[-1]] + DS))) - DS/2
    NT = len(ot)
    NS = len(sedges)

    # calculate Q(s) and Q_s(s)
    Qv=np.zeros((NT, NS))
    Qs=np.zeros((NT, NS))
    # Note that these are organized low s to high s, but still follow
    # the TEF formal definitions from MacCready (2011)
    Qv[:,:-1] = np.fliplr(np.cumsum(np.fliplr(tef_q_lp), axis=1))
    Qs[:,:-1] = np.fliplr(np.cumsum(np.fliplr(tef_qs_lp), axis=1))

    #get bulk values
    Qins=[]
    Qouts=[]
    sins=[]
    souts=[]

    # prepare arrays to hold multi-layer output
    nlay = 30
    QQ = np.nan * np.ones((NT, nlay))
    SS = np.nan * np.ones((NT, nlay))

    if testing:
        plt.close('all')
        dd_list = [154, 260]
        print_info = True
    else:
        dd_list = range(NT)
        print_info = False

    for dd in dd_list:
            
        qv = Qv[dd,:]
        qs = Qs[dd,:]
    
        if print_info == True:
            print('\n**** dd = %d ***' % (dd))
        
        Q_in_m, Q_out_m, s_in_m, s_out_m, div_sal, ind, minmax = tfl.calc_bulk_values(sedges,
            qv, qs, print_info=print_info)
        
        if print_info == True:
            print(' ind = %s' % (str(ind)))
            print(' minmax = %s' % (str(minmax)))
            print(' div_sal = %s' % (str(div_sal)))
            print(' Q_in_m = %s' % (str(Q_in_m)))
            print(' s_in_m = %s' % (str(s_in_m)))
            print(' Q_out_m = %s' % (str(Q_out_m)))
            print(' s_out_m = %s' % (str(s_out_m)))
        
            fig = plt.figure(figsize=(12,8))
        
            ax = fig.add_subplot(121)
            ax.plot(Qv[dd,:], sedges,'.k')
            min_mask = minmax=='min'
            max_mask = minmax=='max'
            print(min_mask)
            print(max_mask)
            ax.plot(Qv[dd,ind[min_mask]], sedges[ind[min_mask]],'*b')
            ax.plot(Qv[dd,ind[max_mask]], sedges[ind[max_mask]],'*r')
            ax.grid(True)
            ax.set_title('Q(s) Time index = %d' % (dd))
            ax.set_ylim(-.1,36.1)
            ax.set_ylabel('Salinity')
        
            ax = fig.add_subplot(122)
            ax.plot(tef_q_lp[dd,:], sbins)
            ax.grid(True)
            ax.set_title('-dQ/ds')
        
        # save multi-layer output
        qq = np.concatenate((Q_in_m, Q_out_m))
        ss = np.concatenate((s_in_m, s_out_m))
        ii = np.argsort(ss)
        if len(ii)>0:
            ss = ss[ii]
            qq = qq[ii]
            NL = len(qq)
            QQ[dd, :NL] = qq
            SS[dd, :NL] = ss

        dd+=1
    
    if testing == False:
        # save results
        bulk = dict()
        bulk['QQ'] = QQ
        bulk['SS'] = SS
        bulk['ot'] = ot
        bulk['qnet_lp'] = qnet_lp
        bulk['fnet_lp'] = fnet_lp
        bulk['ssh_lp'] = ssh_lp
        pickle.dump(bulk, open(out_fn, 'wb'))
    else:
        plt.show()



