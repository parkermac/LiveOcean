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

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'extract/'

# choose the processed TEF file to plot
print('\n%s\n' % '** Choose processed TEF file to plot **')
m_list_raw = os.listdir(indir)
m_list_raw.sort()
m_list = [m for m in m_list_raw if (('.p' in m) and ('tef_' in m))]
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
if False:
    my_npt = int(input('-- Input number -- '))
else:
    my_npt = 1 # for testing
tef_file = m_dict[my_npt]
fn = indir + tef_file

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
if False:
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
plt.close('all')
#fig = plt.figure(figsize=(24,10))

rows = 2
cols = 2

fig, axes = plt.subplots(nrows=rows, ncols=cols, sharex=True, figsize=(24,10))

ax = axes[0,0]
cs = ax.pcolormesh(td, sbins, tef_q_lp.T/1e3, cmap='jet', vmin=-10, vmax=10)
fig.colorbar(cs, ax=ax)
ax.set_xlim(0,365)
ax.set_ylim(35, 25)
ax.set_xlabel('Days')
ax.set_ylabel('Salinity')
ax.text(.05, .9, 'dQ/ds (1e3 m3/s psu-1)', transform=ax.transAxes)
# integrated TEF quantities
#ax.plot(td, Sin, '-k', td, Sout, '--k')
ax.plot(td, S, 'ok', markersize=1)

ax = axes[0,1]
nt, ns = S.shape
Td = np.tile(td.reshape(nt,1),(1, ns))
cs = ax.scatter(Td, Q/1e3, c=S, cmap='jet', vmin=25, vmax=35, s=15)
fig.colorbar(cs, ax=ax)
ax.plot(td, Q[:,0]/1e3, 'g-') # out
ax.plot(td, Q[:,1]/1e3, 'r-') # in
ax.plot(td, Q[:,2]/1e3, 'b-') # out
ax.plot(td, Q[:,3]/1e3, 'm-') # in
ax.set_xlim(0,365)
ax.set_ylim(-400, 400)
ax.set_xlabel('Days')
ax.set_ylabel('Q (1e3 m3/s)')
ax.grid(True)

ax = axes[1,0]
ax.plot(td, qnet_lp/1e3)
ax.set_xlabel('Days')
ax.set_ylim(-50, 50)
ax.set_xlim(0,365)
ax.set_ylabel('LP Volume Flux (1e3 m3/s)')
ax.grid(True)

ax = axes[1,1]
nt, ns = S.shape
Td = np.tile(td.reshape(nt,1),(1, ns))
cs = ax.scatter(Td, S, c=Q/1e3, cmap='jet', vmin=-100, vmax=100, s=15)
fig.colorbar(cs, ax=ax)
ax.plot(td, S[:,0], 'g-') # out
ax.plot(td, S[:,1], 'r-') # in
ax.plot(td, S[:,2], 'b-') # out
ax.plot(td, S[:,3], 'm-') # in
ax.set_xlim(0,365)
ax.set_ylim(35, 25)
ax.set_xlabel('Days')
ax.set_ylabel('S')
ax.grid(True)

# ax = fig.add_subplot(rows, cols, 4)
# ax.plot(td, fnet_lp/1e9)
# ax.set_xlim(0,365)
# ax.set_xlabel('Days')
# ax.set_ylim(-15, 15)
# ax.set_ylabel('Energy Flux (GW)')
# ax.grid(True)

# ax = fig.add_subplot(rows, cols, 5)
# for tt in [10, 50, 100]:
#     # we use these masks because there are multiple values
#     maxmask = Imax[0]==tt
#     minmask = Imin[0]==tt
#     imax = Imax[1][maxmask]
#     imin = Imin[1][minmask]
#     # drop extrema indices which are too close to the ends
#     if len(imax) > 0:
#         mask = np.abs(qcs[tt,imax] - qcs[tt,0]) > crit
#         imax = imax[mask]
#     if len(imax) > 0:
#         mask = np.abs(qcs[tt,imax] - qcs[tt,-1]) > crit
#         imax = imax[mask]
#     if len(imin) > 0:
#         mask = np.abs(qcs[tt,imin] - qcs[tt,0]) > crit
#         imin = imin[mask]
#     if len(imin) > 0:
#         mask = np.abs(qcs[tt,imin] - qcs[tt,-1]) > crit
#         imin = imin[mask]
#     ivec = np.sort(np.concatenate((np.array([0]), imax, imin, np.array([NS]))))
#     ax.plot(qcs[tt,:], sbins)
#     if len(imax) > 0:
#         ax.plot(qcs[tt,imax], sbins[imax], '*r')
#     if len(imax) > 0:
#         ax.plot(qcs[tt,imin], sbins[imin], 'ob')



plt.show()
