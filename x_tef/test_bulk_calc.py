# -*- coding: utf-8 -*-
"""
Code to test the bulk_calc code.
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

ex_list = ['cas6_v1_lo8_2017.01.01_2017.12.31', 'cas6_v2_lo8_2016.12.15_2017.12.31']
dd_offset_list = [0, 17*24]
dd_offset_dict = dict(zip(ex_list, dd_offset_list))

plt.close('all')
fig = plt.figure(figsize=(11,7))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

counter = 0
for ex in ex_list:
    fn = indir0 + ex + '/processed/ai1.p'
    dd_offset = dd_offset_dict[ex]

    print('\nWorking on ' + ex)

    # load the data file
    tef_ex=pickle.load(open(fn, 'rb'))
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
    # tidal averaging
    tef_q_lp = zfun.filt_godin_mat(tef_q)
    tef_qs_lp = zfun.filt_godin_mat(tef_qs)
    qnet_lp = zfun.filt_godin(qnet)
    fnet_lp = zfun.filt_godin(fnet)
    ssh_lp = zfun.filt_godin(ssh)
    pad = 36

    # subsample and cut off nans
    tef_q_lp = tef_q_lp[pad + dd_offset:-(pad+1):24, :]
    tef_qs_lp = tef_qs_lp[pad + dd_offset:-(pad+1):24, :]
    ot = ot[pad + dd_offset:-(pad+1):24]
    qnet_lp = qnet_lp[pad + dd_offset:-(pad+1):24]
    fnet_lp = fnet_lp[pad + dd_offset:-(pad+1):24]
    ssh_lp = ssh_lp[pad + dd_offset:-(pad+1):24]
    
    if counter == 0:
        q0 = tef_q_lp.copy()
    elif counter == 1:
        q1 = tef_q_lp.copy()

    # print(Lfun.modtime_to_datetime(ot[0]))
    # print(Lfun.modtime_to_datetime(ot[-1]))

    # get sizes and make sedges (the edges of sbins)
    DS=sbins[1]-sbins[0]
    sedges = np.concatenate((sbins,np.array([sbins[-1]] + DS))) - DS/2
    NT = len(ot)
    NS = len(sedges)

    # calculate Qv(s) and Qs(s)
    Qv=np.zeros((NT, NS))
    Qs=np.zeros((NT, NS))
    # Note that these are organized low s to high s, but still follow
    # the TEF formal definitions from MacCready (2011)
    #print(tef_q_lp.shape)
    Qv[:,:-1] = np.fliplr(np.cumsum(np.fliplr(tef_q_lp), axis=1))
    Qs[:,:-1] = np.fliplr(np.cumsum(np.fliplr(tef_qs_lp), axis=1))
    
    ax1.plot(Qs[:,0], label=ex)
    ax2.plot(Qv[:,0], label=ex)
    ax3.plot(qnet_lp, label=ex)
    
    counter += 1
    
ax1.legend()
plt.show()


