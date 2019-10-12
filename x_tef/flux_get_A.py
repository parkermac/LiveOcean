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

import flux_fun
reload(flux_fun)

testing = False
verbose = False

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

df = pd.read_pickle(indir + 'flux/two_layer.p')

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
segs = flux_fun.segs

seg_name_list_full = list(segs.keys())
if testing == True:
    seg_name_list = ['J4']
else:
     seg_name_list = seg_name_list_full

# initialize a DataFrame to store results, filled with zeros
df_up_down = pd.DataFrame(0, index=seg_name_list, columns=['q_up', 'q_down', 'A_bad', 'N', 'has_riv'])

# initialize a dict to store output matrices
A_dict = {}

# initialize a DataFrame to store the fluxes used to drive the flux engine
sl2 = []
for sn in seg_name_list:
    sl2.append(sn + '_s')
    sl2.append(sn + '_f')
sl2f = []
for sn in seg_name_list_full:
    sl2f.append(sn + '_s')
    sl2f.append(sn + '_f')
sl2or = ['ocean_s', 'ocean_f'] + ['river_s', 'river_f'] + sl2f
#
q_df = pd.DataFrame(0, index = sl2, columns = sl2or)

for seg_name in seg_name_list:
    
    if verbose:
        print('\n' + seg_name + ':')
    
    seg = segs[seg_name]
    has_riv = False

    # find number of sections for this segment
    N = 0
    for side in list('SNWE'):
        sect_list = seg[side]
        N += len(sect_list) # count up the TEF sections
    if len(seg['R']) > 0:
        N += 1 # the sum of all rivers counts as a single section
        
    # initialize result arrays
    # - inflows
    q = np.nan * np.ones(N)
    f = np.nan * np.ones(N)
    s = np.nan * np.ones(N)
    # - outflows
    Q = np.nan * np.ones(N)
    # Qf is a modified version of Q, non-zero only in the fresher of two layers
    Qf = np.zeros(N)
    F = np.nan * np.ones(N)
    S = np.nan * np.ones(N)
    
    # initialize a list of whether inflow is salt or fresh
    q_sal = list(N*'x')
    # initialize a list of whether outflow is salt or fresh
    Q_sal = list(N*'x')
    
    # initialize a list of section names
    snl = list(N*'x') # the 'x' is just a placeholder
    
    # initialize arrays
    counter = 0
    for side in list('SNWER'):
        sect_list = seg[side]
        if len(sect_list) > 0:
            if side in ['W', 'S']:
                for sect in sect_list:
                    snl[counter] = sect
                    q_s, q_f, f_s, f_f, s_s, s_f, lon, lat = df.loc[sect,:]
                    if q_s > 0:
                        # inflow is salty
                        q[counter] = q_s
                        Q[counter] = -q_f
                        q_sal[counter] = 's'
                        Q_sal[counter] = 'f'
                        Qf[counter] = -q_f
                        f[counter] = f_s
                        F[counter] = -f_f
                        s[counter] = s_s
                        S[counter] = s_f
                    elif q_s < 0:
                        # inflow is fresh
                        q[counter] = q_f
                        Q[counter] = -q_s
                        q_sal[counter] = 'f'
                        Q_sal[counter] = 's'
                        f[counter] = f_f
                        F[counter] = -f_s
                        s[counter] = s_f
                        S[counter] = s_s
                    counter += 1
            elif side in ['E', 'N']:
                for sect in sect_list:
                    snl[counter] = sect
                    q_s, q_f, f_s, f_f, s_s, s_f, lon, lat = df.loc[sect,:]
                    if q_s > 0:
                        # outflow is salty
                        q[counter] =-q_f
                        Q[counter] = q_s
                        q_sal[counter] = 'f'
                        Q_sal[counter] = 's'
                        f[counter] = -f_f
                        F[counter] = f_s
                        s[counter] = s_f
                        S[counter] = s_s
                        
                    elif q_s < 0:
                        # outflow is fresh
                        q[counter] = -q_s
                        Q[counter] = q_f
                        q_sal[counter] = 's'
                        Q_sal[counter] = 'f'
                        Qf[counter] = q_f
                        f[counter] = -f_s
                        F[counter] = f_f
                        s[counter] = s_s
                        S[counter] = s_f
                    counter += 1
            elif side=='R':
                has_riv = True
                Qr = rmean[sect_list].sum()
                q[counter] = Qr
                Q[counter] = 0
                q_sal[counter] = 'f'
                Q_sal[counter] = 'x' # indicate there is no outflow from a river
                f[counter] = 0
                F[counter] = 0
                snl[counter]='river'
    
    if True:
        # make adjustments to enforce volume and salt conservation
        if has_riv:
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
        # then adjust the salinities to match
        s = f/q
        if has_riv:
            S = np.zeros(N)
            S[:-1] = F[:-1]/Q[:-1]
        else:
            S = F/Q
    
        if verbose:    
            print(' - Volume Flux adjustment = %0.5f' % (dq/Q.sum()))
            print(' -   Salt Flux adjustment = %0.5f' % (dqs/F.sum()))

    # The linear system of equations we will solve looks like
    # C * x = y => x = C-1 * y
    # where x is the vector of exchange fractions we are solving for,
    # C is the matrix of incoming fluxes,
    # and y is the vector of outgoing fluxes
    
    if N == 1: # H8 case, only a single section, no rivers
        A = np.array([0], ndmin=2)
        
    else:
        # initialize the system matrices
        C = np.zeros((2*N, 2*N))
        x = np.nan * np.ones(2*N)
        y = np.nan * np.ones(2*N)
        A = np.nan * np.ones((N,N))
        # A is a matrix form of the elements in x, organized as
        # A[incoming section number, outgoing section number], and
        # the vector x is A[:2,:].flatten(order='F')
    
        # fill the elements of C
        Cqf = np.array([ [q[0],q[1]],[f[0],f[1]] ])
        for ii in range(N):
            C[2*ii:2*ii+2, 2*ii:2*ii+2] = Cqf
        
        # make guesses for A
        for cc in range(N):
            A[:N-1,cc] = Q[cc]/Q.sum()
            A[-1,cc] = Qf[cc]/Qf.sum()
                
        def get_A(N, q, Q, f, F, y, A):
            # fill the elements of y:
            for rr in range(N):
                y[2*rr] = Q[rr] - np.sum(A[2:,rr]*q[2:])
                y[2*rr + 1] = F[rr] - np.sum(A[2:,rr]*f[2:])
            # and find the solution
            x = la.tensorsolve(C,y)
            # iCC = la.inv(CC); xalt = iCC.dot(y) # is the same
            # fill in the A matrix with the correct values
            if N > 2:
                A[:2,:] = x.reshape((2,N), order='F')
            else:
                A[:] = x.reshape((N,N), order='F')
            return A
        
        yorig = y.copy()
            
        if (N == 3):# and has_riv:
            # iterative method to scan over all possible river distributions
            # (or section 3 distributions in the case of M4)
            # and then return the average of all valid A's.
            NN = 21 # number of guesses
            AA = np.nan * np.ones((NN,N,N))
            count = 0
            for a20 in np.linspace(0,1,NN):
                A[-1,:] = np.array([a20, 1-a20, 0])
                AA[count,:,:] = get_A(N, q, Q, f, F, yorig, A)
                count += 1
                
            good_guess = []
            for gg in range(NN):
                if (AA[gg,:,:]>=0).all():
                    good_guess.append(gg)
            if verbose:
                print('good guess = ' + str(good_guess))
            if len(good_guess) > 0:
                A = AA[good_guess,:,:].mean(axis=0)
            else: # G5 had no valid solutions
                A = np.ones((N,N)) - np.eye(N)
        else:
            A = get_A(N, q, Q, f, F, yorig, A)
        
        if (N == 2) and (not has_riv):
            # check for bad solutions and force them to have no mixing
            if (A < 0).any():
                if verbose:
                    print('Bad value for N=2 segment ' + seg_name)
                    print('ORIGINAL')
                    np.set_printoptions(precision=3, suppress=True)
                    print(A)
                    print('NEW')
                A = np.ones((N,N)) - np.eye(N)
                if verbose:
                    np.set_printoptions(precision=3, suppress=True)
                    print(A)
                    
        if (N == 2) and has_riv:
            A = np.ones((N,N)) - np.eye(N)
            
    # store results in the DataFrame
    for ii in range(N-1):
        this_q = q[ii] * A[ii,ii]
        this_up = df_up_down.loc[seg_name, 'q_up']
        this_down = df_up_down.loc[seg_name, 'q_down']
        if (q_sal[ii] == 's') and (this_q > this_up):
            df_up_down.loc[seg_name, 'q_up'] = this_q
        if (q_sal[ii] == 'f') and (this_q > this_down):
            df_up_down.loc[seg_name, 'q_down'] = this_q
    if (A < 0).any():
        df_up_down.loc[seg_name, 'A_bad'] = 1
    df_up_down.loc[seg_name, 'N'] = N
    if has_riv:
        df_up_down.loc[seg_name, 'has_riv'] = 1
        
    # store other results in dicts
    A_dict[seg_name] = A
    
    # *************************************************
    # store results in the DataFrame of fluxes used for the flux engine
    def update_q(q_df, q_sal, seg_name, Seg_name, q, Q, sss):
        # transport in from adjoining segment
        # and vertical transport up from below or down from above
        if q_sal[sss] == 's':
            q_df.loc[seg_name+'_s',Seg_name+'_s'] += q[sss]*(1-A[sss,sss])
            q_df.loc[seg_name+'_f',Seg_name+'_s'] += q[sss] * A[sss,sss]
        elif q_sal[sss] == 'f':
            q_df.loc[seg_name+'_f',Seg_name+'_f'] += q[sss]*(1-A[sss,sss])
            q_df.loc[seg_name+'_s',Seg_name+'_f'] += q[sss] * A[sss,sss]
        # transport out to adjoining segment
        # in this case the column (where tracer comes from) is the segment itself
        if Q_sal[sss] == 's':
            q_df.loc[seg_name+'_s',seg_name+'_s'] -=  Q[sss]
        elif Q_sal[sss] == 'f':
            q_df.loc[seg_name+'_f',seg_name+'_f'] -= Q[sss]
        return q_df
        
    sss = 0
    for sect in snl:
        if (seg_name=='J1' and sect=='jdf1') or (seg_name=='G6' and sect=='sog5'):
            Seg_name = 'ocean'
            if verbose:
                print(' -- Adjoining segment = ' + Seg_name + ' for ' + sect)
            q_df = update_q(q_df, q_sal, seg_name, Seg_name, q, Q, sss)
        elif sect == 'river':
            q_df.loc[seg_name+'_f','river_f'] = q[sss]
        else:
            # find which segment each section connects to
            for Seg_name in segs.keys():
                Seg = segs[Seg_name]
                for side in list('SNWE'):
                    sect_list = Seg[side]
                    if sect in sect_list and Seg_name != seg_name:
                        if verbose:
                            print(' -- Adjoining segment = ' + Seg_name + ' for ' + sect)
                        q_df = update_q(q_df, q_sal, seg_name, Seg_name, q, Q, sss)
        sss += 1
        
    # add net vertical advection in the section (required for W4)
    qx_s = q_df.loc[seg_name+'_s'].sum(axis=0)
    qx_f = q_df.loc[seg_name+'_f'].sum(axis=0)
    if np.abs(qx_s + qx_f) > 1:
        print('qx warning')
    if qx_s > 0:
        # upward transport
        q_df.loc[seg_name+'_s',seg_name+'_s'] -=  qx_s
        q_df.loc[seg_name+'_f',seg_name+'_s'] +=  qx_s
    elif qx_s <= 0:
        # downward transport
        q_df.loc[seg_name+'_f',seg_name+'_f'] +=  qx_s
        q_df.loc[seg_name+'_s',seg_name+'_f'] -=  qx_s
    
    if verbose:
        print('Solution matrix A:')
        np.set_printoptions(precision=3, suppress=True)
        print(A)
        print(q_sal)
        print('')
        
    if not testing:
        # save results
        q_df.to_pickle(indir + 'flux/q_df.p')
        

