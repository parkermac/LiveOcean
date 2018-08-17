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

# choose the processed TEF fiel to plot
print('\n%s\n' % '** Choose processed TEF file to process **')
m_list_raw = os.listdir(indir)
m_list_raw.sort()
m_list = [m for m in m_list_raw if (('.p' in m) and ('tef_' in m))]
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
if True:
    my_npt = int(input('-- Input number -- '))
else:
    my_npt = 0 # for testing
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
tef_q = zfun.filt_godin_mat(tef_q)
tef_qs = zfun.filt_godin_mat(tef_qs)
qnet_lp = zfun.filt_godin(qnet)
fnet_lp = zfun.filt_godin(fnet)

# subsample
tef_q = tef_q[36:-37:24, :]
tef_qs = tef_qs[36:-37:24, :]
td = td[36:-37:24]
qnet_lp = qnet_lp[36:-37:24]
fnet_lp = fnet_lp[36:-37:24]

# # find integrated TEF quantities
# # alternate method using cumulative sum of the transport
# # to identify the salinity dividing inflow and outflow
# # RESULT: this way is not sensitive to the number of
# # salinity bins.
# #
# start by making the low-passed flux arrays sorted
# from high to low salinity
rq = np.fliplr(tef_q)
rqs = np.fliplr(tef_qs)
# then form the cumulative sum (the function Q(s))
qcs = np.cumsum(rq, axis=1)
nt = len(td)


# new version to handle more layers
from scipy.signal import argrelextrema
Imax = argrelextrema(qcs, np.greater, axis=1, order=int(NS/50))
Imin = argrelextrema(qcs, np.less, axis=1, order=int(NS/50))

nlay = 10
Q = np.zeros((nt, nlay))
QS = np.zeros((nt, nlay))

for tt in range(nt):
    maxmask = Imax[0]==tt
    minmask = Imin[0]==tt
    imax = Imax[1][maxmask]
    imin = Imin[1][minmask]
    ivec = (np.concatenate((np.array([0]), imax, imin, np.array([NS])))) # use np.sort()?
    for ii in range(len(ivec)-1):
        Q[tt,ii] = rq[tt, ivec[ii]:ivec[ii+1]].sum()
        QS[tt,ii] = rqs[tt, ivec[ii]:ivec[ii+1]].sum()

# form derived quantities
Q[Q==0] = np.nan
S = QS/Q

# plotting
plt.close('all')
fig = plt.figure(figsize=(16,10))

ax = fig.add_subplot(221)
cs = ax.pcolormesh(td, sbins, tef_q.T/1e3, cmap='jet', vmin=-10, vmax=10)
fig.colorbar(cs, ax=ax)
ax.set_xlim(0,365)
ax.set_ylim(35, 20)
ax.set_xlabel('Days')
ax.set_ylabel('Salinity')
ax.text(.05, .9, 'dQ/ds (1e3 m3/s psu-1)', transform=ax.transAxes)
# integrated TEF quantities
#ax.plot(td, Sin, '-k', td, Sout, '--k')
ax.plot(td, S, '*k')

ax = fig.add_subplot(222)
ax.plot(td, Q/1e3,)
ax.set_xlim(0,365)
ax.set_ylim(-500, 500)
ax.set_xlabel('Days')
ax.set_ylabel('Q (1e3 m3/s)')


ax = fig.add_subplot(223)
ax.plot(td, qnet_lp/1e3)
ax.set_xlabel('Days')
ax.set_ylim(-20, 50)
ax.set_xlim(0,365)
ax.set_ylabel('LP Volume Flux (1e3 m3/s)')

ax = fig.add_subplot(224)
ax.plot(td, fnet_lp/1e9)
ax.set_xlim(0,365)
ax.set_xlabel('Days')
ax.set_ylim(0, 15)
ax.set_ylabel('Energy Flux (GW)')


plt.show()
