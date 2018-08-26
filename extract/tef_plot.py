"""
Process a TEF extraction.
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

import tef_fun
from importlib import reload
reload(tef_fun)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()


from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

Ldir = Lfun.Lstart()

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

indir = Ldir['LOo'] + 'tef/'
print('\nSelect an Extraction to plot:\n')
List_raw = os.listdir(indir)
List_raw.sort()
List = [item for item in List_raw]
NL = len(List)
Ldict = dict(zip(range(NL), List))
for ii in range(NL):
    print(str(ii) + ': ' + List[ii])
if True:
    my_ii = int(input('-- Input number: '))
else:
    my_ii = 0 # for testing
Litem = Ldict[my_ii]
print('\nProcessing ' + Litem + '\n')
LList_raw = os.listdir(indir + Litem)
LList_raw.sort()
LList = [item for item in LList_raw if ('.p' in item)]
Indir = indir + Litem + '/'

for LL in ['AInorth.p','AIsouth.p', 'MBnorth.p', 'MBmidnorth.p', 'MBmid.p']: #LList:
    
    sect_name = LL.replace('.p','')
    print('\n** ' + sect_name + ' **')

    fn = Indir + LL

    # load results
    tef_dict = pickle.load(open(fn, 'rb'))
    tef_q = tef_dict['tef_q']
    tef_qs = tef_dict['tef_qs']
    sbins = tef_dict['sbins']
    qnet = tef_dict['qnet']
    fnet = tef_dict['fnet']
    ot = tef_dict['ot']
    td = (ot - ot[0])/86400
    NS = len(sbins)

    # low-pass
    if True:
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

    # # find integrated TEF quantities
    # # alternate method using cumulative sum of the transport
    # # to identify the salinity dividing inflow and outflow
    # # RESULT: this way is not sensitive to the number of
    # # salinity bins.
    # #
    # start by making the low-passed flux arrays sorted
    # from high to low salinity
    rq = np.fliplr(tef_q_lp)
    rqs = np.fliplr(tef_qs_lp)
    # then form the cumulative sum (the function Q(s))
    qcs = np.cumsum(rq, axis=1)
    nt = len(td)

    # new version to handle more layers
    from scipy.signal import argrelextrema
    Imax = argrelextrema(qcs, np.greater, axis=1, order=int(NS/50))
    Imin = argrelextrema(qcs, np.less, axis=1, order=int(NS/50))

    nlay_max = 4
    Q = np.zeros((nt, nlay_max))
    QS = np.zeros((nt, nlay_max))

    crit = np.nanmax(np.abs(qcs)) / 50

    for tt in range(nt):
        # we use these masks because there are multiple values for a given day
        maxmask = Imax[0]==tt
        minmask = Imin[0]==tt
        imax = Imax[1][maxmask]
        imin = Imin[1][minmask]
        # drop extrema indices which are too close to the ends
        if len(imax) > 0:
            mask = np.abs(qcs[tt,imax] - qcs[tt,0]) > crit
            imax = imax[mask]
        if len(imax) > 0:
            mask = np.abs(qcs[tt,imax] - qcs[tt,-1]) > crit
            imax = imax[mask]
        if len(imin) > 0:
            mask = np.abs(qcs[tt,imin] - qcs[tt,0]) > crit
            imin = imin[mask]
        if len(imin) > 0:
            mask = np.abs(qcs[tt,imin] - qcs[tt,-1]) > crit
            imin = imin[mask]
        ivec = np.sort(np.concatenate((np.array([0]), imax, imin, np.array([NS]))))
        nlay = len(ivec)-1
    
        # combine non-alternating layers
        qq = np.zeros(nlay)
        qqs = np.zeros(nlay)
        jj = 0
        for ii in range(nlay):
            qlay = rq[tt, ivec[ii]:ivec[ii+1]].sum()
            qslay = rqs[tt, ivec[ii]:ivec[ii+1]].sum()
            if ii == 0:
                qq[0] = qlay
                qqs[0] = qslay
            else:
                if np.sign(qlay)==np.sign(qq[jj]):
                    qq[jj] += qlay
                    qqs[jj] += qslay
                    nlay -= 1
                else:
                    jj += 1
                    qq[jj] = qlay
                    qqs[jj] = qslay

        if nlay == 1:
            if qq[0] >= 0:
                Q[tt,1] = qq[0]
                QS[tt,1] = qqs[0]
            elif qq[0] < 0:
                Q[tt,2] = qq[0]
                QS[tt,2] = qqs[0]
        elif nlay ==2:
            if qq[0] >= 0:
                Q[tt,1] = qq[0]
                QS[tt,1] = qqs[0]
                Q[tt,2] = qq[1]
                QS[tt,2] = qqs[1]
            elif qq[0] < 0:
                Q[tt,1] = qq[1]
                QS[tt,1] = qqs[1]
                Q[tt,2] = qq[0]
                QS[tt,2] = qqs[0]
        elif nlay ==3:
            if qq[0] >= 0:
                Q[tt,1] = qq[0]
                QS[tt,1] = qqs[0]
                Q[tt,2] = qq[1]
                QS[tt,2] = qqs[1]
                Q[tt,3] = qq[2]
                QS[tt,3] = qqs[2]
            elif qq[0] < 0:
                Q[tt,0] = qq[0]
                QS[tt,0] = qqs[0]
                Q[tt,1] = qq[1]
                QS[tt,1] = qqs[1]
                Q[tt,2] = qq[2]
                QS[tt,2] = qqs[2]
        elif nlay ==4:
            if qq[0] >= 0:
                print('- Backwards 4: td=%5.1f  nlay = %d' % (td[tt],nlay))
                Q[tt,:] = np.nan
                QS[tt,:] = np.nan
            elif qq[0] < 0:
                Q[tt,0] = qq[0]
                QS[tt,0] = qqs[0]
                Q[tt,1] = qq[1]
                QS[tt,1] = qqs[1]
                Q[tt,2] = qq[2]
                QS[tt,2] = qqs[2]
                Q[tt,3] = qq[3]
                QS[tt,3] = qqs[3]
        else:
            print('- Excess layers: td=%5.1f  nlay = %d' % (td[tt],nlay))
            Q[tt,:] = np.nan
            QS[tt,:] = np.nan
        


    # form derived quantities
    Q[Q==0] = np.nan
    S = QS/Q

    # plotting
    #plt.close('all')
    #fig = plt.figure(figsize=(24,10))

    lw = 2.5

    fig, axes = plt.subplots(nrows=2, ncols=3, sharex=True, figsize=(20,8))

    ax = axes[0,0]
    nt, ns = S.shape
    Td = np.tile(td.reshape(nt,1),(1, ns))
    ax.plot(td, Q[:,3]/1e3, '-r', linewidth=lw, label='Qin shallow', alpha=.5)
    ax.plot(td, Q[:,2]/1e3, '-b', linewidth=lw, label='Qout') # out
    ax.plot(td, Q[:,1]/1e3, '-r', linewidth=lw, label='Qin') # in
    ax.plot(td, Q[:,0]/1e3, '-b', linewidth=lw, label='Qout deep', alpha=.5) # out
    ax.legend(ncol=2, loc='upper left')
    ax.set_xlim(0,365)
    # ax.set_ylim(-500, 500)
    ax.set_ylim(-50, 50)
    ax.set_xlabel('Days')
    ax.set_ylabel('Q (1e3 m3/s)')
    ax.grid(True)
    
    print('^^ %s Qin_mean = %0.1f 1e3 m3/s' % (sect_name, np.nanmean(Q[:,1]/1e3)))
    print('^^ %s QSin_mean = %0.1f psu 1e3 m3/s' % (sect_name, np.nanmean(QS[:,1]/1e3)))

    ax = axes[1,0]
    nt, ns = S.shape
    Td = np.tile(td.reshape(nt,1),(1, ns))
    ax.plot(td, S[:,3], '-r', linewidth=lw, label='Sin shallow', alpha=.5)
    ax.plot(td, S[:,2], '-b', linewidth=lw, label='Sout') # out
    ax.plot(td, S[:,1], '-r', linewidth=lw, label='Sin') # in
    ax.plot(td, S[:,0], '-b', linewidth=lw, label='Sout deep', alpha=.5) # out
    # mark reversals
    rev = (S[:,2] - S[:,1]) > 0
    ax.plot(td[rev], 34.8*np.ones(nt)[rev], '*k', label='Reversals')
    ax.legend(loc='upper right')
    ax.set_xlim(0,365)
    ax.set_ylim(35, 10)
    ax.set_xlabel('Days')
    ax.set_ylabel('Salinity')
    ax.grid(True)

    ax = axes[0,1]
    ax.plot(td, qnet_lp/1e3, '-k', linewidth=lw)
    ax.set_xlabel('Days')
    # ax.set_ylim(-50, 50)
    ax.set_ylim(-5, 5)
    ax.set_xlim(0,365)
    ax.set_ylabel('LP Volume Flux (1e3 m3/s)')
    ax.grid(True)
    ax.text(0.05, 0.1, 'Positive is Landward', transform=ax.transAxes)

    ax = axes[1,1]
    ax.plot(td, fnet_lp/1e9, '-k', linewidth=lw)
    ax.set_xlim(0,365)
    ax.set_xlabel('Days')
    #ax.set_ylim(-15, 15)
    ax.set_ylim(-1, 1)
    ax.set_ylabel('Energy Flux (GW)')
    ax.grid(True)

    # remove extra axes
    axes[0,2].set_axis_off()
    axes[1,2].set_axis_off()

    # add section location map
    ax = fig.add_subplot(1,3,3)
    
    x0, x1, y0, y1, landward = sect_df.loc[sect_name,:]
    ax.plot([x0, x1], [y0, y1], '-m', linewidth=3)
    ax.set_title(sect_name)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.set_xlim(-124, -122)
    ax.set_ylim(47, 49)


    plt.show()
