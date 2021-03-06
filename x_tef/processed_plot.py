"""
Plot processed fluxes as a pcolor (time, salinity), with the goal of learning
about mixing and transport at sills.

"""
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
Ldir = Lfun.Lstart()

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

import numpy as np
import matplotlib.pyplot as plt
import pickle
from datetime import datetime, timedelta

import tef_fun
# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

# from warnings import filterwarnings
# filterwarnings('ignore') # skip some warning messagesd

# ===================================================================
testing = True
year = 2017
year_str = str(year)
# ===================================================================

# choose input and organize output
Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'tef/'
indir0 = indir0 + 'cas6_v3_lo8b_'+year_str+'.01.01_'+year_str+'.12.31' + '/'
indir = indir0 + 'processed/'

sect_list = list(sect_df.index)
if testing:
    sect_list = ['ai2','tn2']

plt.close('all')

for sn in sect_list:
    bulk = pickle.load(open(indir + sn + '.p', 'rb'))
    # bulk is a dict with keys ['tef_q', 'tef_qs', 'sbins', 'ot', 'qnet', 'fnet', 'ssh']
    ot = bulk['ot'] # model time in seconds from 1/1/1970 (hourly)
    sbins = bulk['sbins']
    q = bulk['tef_q'] # packed (time, salinity bin)
    qf = zfun.filt_godin_mat(q)
    qnet = bulk['qnet']
    day = (ot - ot[0])/86400
    
    # plotting
    fig = plt.figure(figsize=(14,10))
    t0 = 75
    t1 = 95
    slo = 26
    shi = 34
    
    islo = (np.abs(sbins-slo)).argmin()
    ishi = (np.abs(sbins-shi)).argmin()
    it0 = (np.abs(day - t0)).argmin()
    it1 = (np.abs(day - t1)).argmin()
    sbins = sbins[islo:ishi]
    day = day[it0:it1]
    qnet = qnet[it0:it1]
    q = q[it0:it1, islo:ishi]
    qf = qf[it0:it1, islo:ishi]
    
    ax = fig.add_subplot(311)
    scl = qnet.std()/20
    cs = ax.pcolormesh(day, sbins, q.T, cmap='bwr', vmin=-scl, vmax=scl)
    #fig.colorbar(cs, ax=ax)
    ax.axis([t0, t1, slo, shi])
    ax.grid(True)
    ax.set_ylabel('Salinity')
    ax.set_title(sn)
    ax.text(.05,.8,
    'Raw Transport in Salinity Bins\nColor scale is +/- %d (m3/s / (.001psu))' % (scl),
    transform = ax.transAxes)

    ax = fig.add_subplot(312)
    scl = qnet.std()/100
    cs = ax.pcolormesh(day, sbins, qf.T, cmap='bwr', vmin=-scl, vmax=scl)
    #fig.colorbar(cs, ax=ax)
    ax.axis([t0, t1, slo, shi])
    ax.grid(True)
    ax.set_ylabel('Salinity')
    ax.text(.05,.8,
    'Tidally averaged Transport in Salinity Bins\nColor scale is +/- %d (m3/s / (.001psu))' % (scl),
    transform = ax.transAxes)
    
    ax = fig.add_subplot(313)
    ax.plot(day, qnet/1000, '-', color='g')
    ax.set_xlim(t0, t1)
    ax.grid(True)
    ax.set_xlabel('Yearday')
    ax.set_ylabel('Qnet (1000 m3/s)')
    
plt.show()