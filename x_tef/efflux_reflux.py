"""
Calculate efflux-reflux coefficients.
"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import pickle
import pandas as pd
from datetime import datetime, timedelta

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')


pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import tef_fun
from importlib import reload
reload(tef_fun)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

# select input file
if False:
    indir0 = Ldir['LOo'] + 'tef/'
    # choose the tef extraction to plot
    item = Lfun.choose_item(indir0)
    indir = indir0 + item + '/'
else:
    indir = '/Users/pm7/Documents/LiveOcean_output/tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'

df = pd.read_pickle(indir + 'two_layer/TwoLayer.p')

# get river flow from the model run
fnr = Ldir['gtag'] + '_2017.01.01_2018.12.31.p'
fn = Ldir['LOo'] + 'river/' + fnr
dfr = pd.read_pickle(fn)
dt0 = datetime(2017,1,1)
dt1 = datetime(2017,12,31)
dfrr = dfr.loc[dt0:dt1,:]
rmean = dfrr.mean() # a series


# segment definitions, assembled by looking at the figure
# created by plot_thalweg_mean.py
segs = {
        'J1':{'S':[], 'N':[], 'W':['jdf1'], 'E':['jdf2'], 'R':['sanjuan', 'hoko']},
        #'J2':{'S':[], 'N':[], 'W':['jdf2'], 'E':['jdf3'], 'R':[]},
        }

for seg_name in segs.keys():
    
    seg = segs[seg_name]

    # find number of sections for this segment
    N = 0
    for side in list('NSEW'):
        s = seg[side]
        N += len(s)
    if len(seg['R']) > 0:
        N += 1
    # initialize result arrays
    # - inflows
    q = np.nan * np.ones(N)
    f = np.nan * np.ones(N)
    # - outflows
    Q = np.nan * np.ones(N)
    F = np.nan * np.ones(N)

    # initialize arrays
    counter = 0
    for side in list('WESNR'):
        s = seg[side]
        if len(s) > 0:
            if side in ['W', 'S']:
                for sect in s:
                    q_s, q_f, f_s, f_f, s_s, s_f = df.loc[sect,:]
                    if q_s > 0:
                        q[counter] = q_s
                        Q[counter] = -q_f
                        f[counter] = f_s
                        F[counter] = -f_f
                    elif q_s < 0:
                        q[counter] = q_f
                        Q[counter] = -q_s
                        f[counter] = f_f
                        F[counter] = -f_s
                    counter += 1
            elif side in ['E', 'N']:
                for sect in s:
                    q_s, q_f, f_s, f_f, s_s, s_f = df.loc[sect,:]
                    if q_s > 0:
                        q[counter] =-q_f
                        Q[counter] = q_s
                        f[counter] = -f_f
                        F[counter] = f_s
                    elif q_s < 0:
                        q[counter] = -q_s
                        Q[counter] = q_f
                        f[counter] = -f_s
                        F[counter] = f_f
                    counter += 1
            elif side=='R':
                Qr = rmean[s].sum()
                q[counter] = Qr
                Q[counter] = 0
                f[counter] = 0
                F[counter] = 0
                
    # make adjustments to enforce volume and salt conservation
    if len(seg['R'])>0:
        # if there are rivers only adjust the non-river fluxes
        dq = q.sum() - Q.sum()
        q[:-1] -= 0.5*dq/(N-1)
        Q[:-1] += 0.5*dq/(N-1)
        dqs = f.sum() - F.sum()
        f[:-1] -= 0.5*dqs/(N-1)
        F[:-1] += 0.5*dqs/(N-1)
    else:
        # if there are no rivers adjust all fluxes
        dq = q.sum() - Q.sum()
        q[:] -= 0.5*dq/N
        Q[:] += 0.5*dq/N
        dqs = f.sum() - F.sum()
        f[:] -= 0.5*dqs/N
        F[:] += 0.5*dqs/N
    print('\n' + seg_name + ':')
    print(' - Volume Flux adjustment = %0.5f' % (dq/Q.sum()))
    print(' -   Salt Flux adjustment = %0.5f' % (dqs/F.sum()))

    # The linear system of equations we will solve looks like
    # C * x = y => x = C-1 * y
    # where x is the vector of exchange fractions we are solving for,
    # C is the matrix of incoming fluxes,
    # and y is the vector of outgoing fluxes
    
    # initialize the system matrices
    C = np.zeros((2*N, N*N))
    x = np.nan * np.ones(N*N)
    y = np.nan * np.ones(2*N)
    A = np.nan * np.ones((N,N))
    # A is a matrix form of the elements in x, organized as
    # A[incoming section number, outgoing section number], and
    # the vector x is A.flatten(order='F')
    
    # fill the elements of C
    for rr in range(N):
        for cc in range(N):
            if rr == cc:
                C[rr*2,cc*N:(cc*N+N)] = q
                C[rr*2 + 1,cc*N:(cc*N+N)] = f
                
    # fill the elements of y:
    for rr in range(N):
        y[2*rr] = Q[rr]
        y[2*rr + 1] = F[rr]
    
    # to solve the system we need to make guesses for elements of A,
    # then move those to the RHS and remove columns from C
    
    # The code below may only work for dN = 3
    
    dN = N*N - 2*N # number of guesses needed
    
    a20_vec = np.linspace(.9,1,11)
    for a20 in a20_vec:
    
        if dN > 0:
        
            # make guess for a row of A
        
            # BUT it may be reasonable to assume that the river flow
            # only goes to the upper outgoing layers, not the incoming ones...
        
            for ii in range(dN):
                A[N-1, ii] = Q[ii]/Q.sum()
            
            # Experimentation
            if seg_name == 'J1':
                A[N-1, :] = np.array([a20, 1-a20, 0])
    
            # make adjustments to y
            yy = y.copy()
            for ii in range(dN):
                yy[2*ii] -= A[N-1,ii] * q[N-1]
                yy[2*ii + 1] -= A[N-1,ii] * f[N-1]
        
            # CC is C with appropriate columns removed
            CC = np.nan * np.ones((2*N,2*N))
            for ii in range(dN):
                CC[:,(N-1)*ii] = C[:,N*ii]
                CC[:,(N-1)*ii + 1] = C[:,N*ii + 1]
        else:
            CC = C.copy()
            yy = y.copy()

        # and find the solution
        x = la.tensorsolve(CC,yy)
        # iCC = la.inv(CC)
        # xalt = iCC.dot(y)
    
        # fill in rest of the A matrix
        if dN > 0:
            A[:-1] = x.reshape((N-1,N), order='F')
        else:
            A[:] = x.reshape((N,N), order='F')
        
        print('')
        print(A)
    

