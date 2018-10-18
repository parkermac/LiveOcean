"""
Process a TEF extraction using the new multi-layer
code from Martin (?) Lorenz.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import pickle

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun
Ldir = Lfun.Lstart()

import tef_fun
from importlib import reload
reload(tef_fun)
reload(zfun)

import tef_fun_lorenz as tfl
reload(tfl)

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

indir0 = Ldir['LOo'] + 'tef/'
# choose the tef extraction to plot
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/'

LList_raw = os.listdir(indir )
LList_raw.sort()
LList = [item for item in LList_raw if ('.p' in item) and ('ThalMean' not in item)]

if False: # plot all .p files
    save_fig = True
    out_dir = indir + 'plots/'
    Lfun.make_dir(out_dir)
else: # override
    save_fig = False
    LList = [item for item in LList if 'jdf1' in item]

plt.close('all')

for LL in LList:
    
    sect_name = LL.replace('.p','')
    print('\n** ' + sect_name + ' **')

    fn = indir + LL
    
    tef_dict = pickle.load(open(fn, 'rb'))
    tef_q = tef_dict['tef_q']
    tef_qs = tef_dict['tef_qs']
    sbins = tef_dict['sbins']
    smax = sbins.max()
    qnet = tef_dict['qnet']
    fnet = tef_dict['fnet']
    ot = tef_dict['ot']
    td = (ot - ot[0])/86400
    NS = len(sbins)

    # low-pass
    tidal_average = False
    if tidal_average:
        # tidal averaging
        tef_q_lp = zfun.filt_godin_mat(tef_q)
        tef_qs_lp = zfun.filt_godin_mat(tef_qs)
        qnet_lp = zfun.filt_godin(qnet)
        fnet_lp = zfun.filt_godin(fnet)
        pad = 36
    else:
        # nday Hanning window
        nday = 5
        nfilt = nday*24
        tef_q_lp = zfun.filt_hanning_mat(tef_q, n=nfilt)
        tef_qs_lp = zfun.filt_hanning_mat(tef_qs, n=nfilt)
        qnet_lp = zfun.filt_hanning(qnet, n=nfilt)
        fnet_lp = zfun.filt_hanning(fnet, n=nfilt)
        pad = int(np.ceil(nfilt/2))

    # subsample
    tef_q_lp = tef_q_lp[pad:-(pad+1):24, :]
    tef_qs_lp = tef_qs_lp[pad:-(pad+1):24, :]
    td = td[pad:-(pad+1):24]
    qnet_lp = qnet_lp[pad:-(pad+1):24]
    fnet_lp = fnet_lp[pad:-(pad+1):24]
    
    # form integrals over s (low to high)
    rq = np.fliplr(tef_q_lp)
    rqs = np.fliplr(tef_qs_lp)
    # then form the cumulative sum (the function Q(s))
    Q = np.cumsum(rq, axis=1)
    QS = np.cumsum(rqs, axis=1)
    Qr = np.fliplr(Q)
    QSr = np.fliplr(QS)
    s = sbins
    
    NT = tef_q_lp.shape[0]
    NLmax = 10
    Qin = 0 * np.ones((NT, NLmax))
    Qout = 0 * np.ones((NT, NLmax))
    Sin = 0 * np.ones((NT, NLmax))
    Sout = 0 * np.ones((NT, NLmax))
    for tt in range(NT):
        s = sbins
        Qv = Qr[tt,:]
        Qs = QSr[tt,:]
        comp = 20
        min_trans = 100
    
        Q_in_m,Q_out_m,s_in_m,s_out_m,div_sal = tfl.calc_bulk_values(s,Qv,Qs,comp,min_trans)
        Qin[tt,0:len(Q_in_m)] = Q_in_m
        Qout[tt,0:len(Q_out_m)] = Q_out_m
        Sin[tt,0:len(s_in_m)] = s_in_m
        Sout[tt,0:len(s_out_m)] = s_out_m
        
    plt.close('all')
    fig = plt.figure()
    
    ax = fig.add_subplot(211)
    for ll in range(NLmax):
        ax.plot(Qin[:,ll]/1000)
        ax.plot(Qout[:,ll]/1000)
    
    ax = fig.add_subplot(212)
    for ll in range(NLmax):
        ax.plot(Sin[:,ll])
        ax.plot(Sout[:,ll])
        
    plt.show()
    


