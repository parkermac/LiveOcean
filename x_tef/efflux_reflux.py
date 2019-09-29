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

testing = False

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
        'J2':{'S':[], 'N':[], 'W':['jdf2'], 'E':['jdf3'], 'R':[]},
        'J3':{'S':[], 'N':[], 'W':['jdf3'], 'E':['jdf4'], 'R':['elwha']},
        'J4':{'S':[], 'N':['sji1'], 'W':['jdf4'], 'E':['ai1','dp'], 'R':['dungeness']},
        
        'G1':{'S':['sji1'], 'N':['sji2'], 'W':[], 'E':[], 'R':['samish']},
        'G2':{'S':['sji2'], 'N':['sog1'], 'W':[], 'E':[], 'R':['nooksack', 'cowichan']},
        'G3':{'S':['sog1'], 'N':['sog2'], 'W':[], 'E':[], 'R':['nanaimo', 'fraser']},
        'G4':{'S':['sog2'], 'N':[], 'W':['sog3'], 'E':[], 'R':['clowhom', 'squamish']},
        'G5':{'S':[], 'N':['sog4'], 'W':[], 'E':['sog3'], 'R':['englishman', 'tsolum', 'oyster']},
        'G6':{'S':['sog4'], 'N':['sog5'], 'W':[], 'E':[], 'R':[]},
        
        'M1':{'S':['ai2'], 'N':[], 'W':['ai1'], 'E':[], 'R':[]},
        'M2':{'S':['ai3'], 'N':['ai2'], 'W':[], 'E':[], 'R':[]},
        'M3':{'S':['hc1'], 'N':['ai3'], 'W':[], 'E':['ai4'], 'R':[]},
        'M4':{'S':['mb1'], 'N':['wb1'], 'W':['ai4'], 'E':[], 'R':[]},
        'M5':{'S':['mb2'], 'N':['mb1'], 'W':[], 'E':[], 'R':[]},
        'M6':{'S':['mb3'], 'N':['mb2'], 'W':[], 'E':[], 'R':['green', 'cedar']},
        'M7':{'S':['mb4'], 'N':['mb3'], 'W':[], 'E':[], 'R':[]},
        'M8':{'S':['mb5'], 'N':['mb4'], 'W':[], 'E':[], 'R':[]},
        'M9':{'S':['tn1'], 'N':['mb5'], 'W':[], 'E':[], 'R':['puyallup']},
        
        'S1':{'S':['tn2'], 'N':['tn1'], 'W':[], 'E':[], 'R':[]},
        'S2':{'S':['tn3'], 'N':['tn2'], 'W':[], 'E':[], 'R':[]},
        'S3':{'S':[], 'N':['tn3'], 'W':['ss1'], 'E':[], 'R':[]},
        'S4':{'S':[], 'N':[], 'W':['ss2'], 'E':['ss1'], 'R':['nisqually']},
        'S5':{'S':[], 'N':[], 'W':['ss3'], 'E':['ss2'], 'R':[]},
        'S6':{'S':[], 'N':[], 'W':[], 'E':['ss3'], 'R':['deschutes']},
        
        'W1':{'S':['wb1'], 'N':['wb2'], 'W':[], 'E':[], 'R':['snohomish']},
        'W2':{'S':['wb2'], 'N':['wb3'], 'W':[], 'E':[], 'R':['stillaguamish']},
        'W3':{'S':['wb3'], 'N':[], 'W':[], 'E':['wb4'], 'R':[]},
        'W4':{'S':[], 'N':[], 'W':['wb4', 'dp'], 'E':[], 'R':['skagit']},
        
        'H1':{'S':['hc2'], 'N':['hc1'], 'W':[], 'E':[], 'R':[]},
        'H2':{'S':[], 'N':['hc2'], 'W':['hc3'], 'E':[], 'R':[]},
        'H3':{'S':['hc4'], 'N':[], 'W':[], 'E':['hc3'], 'R':['duckabush', 'dosewallips']},
        'H4':{'S':['hc5'], 'N':['hc4'], 'W':[], 'E':[], 'R':['hamma']},
        'H5':{'S':['hc6'], 'N':['hc5'], 'W':[], 'E':[], 'R':[]},
        'H6':{'S':[], 'N':['hc6'], 'W':[], 'E':['hc7'], 'R':['skokomish']},
        'H7':{'S':[], 'N':[], 'W':['hc7'], 'E':['hc8'], 'R':[]},
        'H8':{'S':[], 'N':[], 'W':['hc8'], 'E':[], 'R':[]},
        
        #'##':{'S':[], 'N':[], 'W':[], 'E':[], 'R':[]},
        }

if testing == True:
    seg_name_list = ['M4']
else:
     seg_name_list = segs.keys()

# initialize a DataFrame to store results, filled with zeros
df_up_down = pd.DataFrame(0, index=seg_name_list, columns=['q_up', 'q_down', 'A_bad', 'N', 'has_riv'])

for seg_name in seg_name_list:
    
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
    
    # initialize a list of whether in inflow is salt or fresh
    q_sal = list(N*'f') 
    
    # initialize arrays
    counter = 0
    for side in list('SNWER'):
    #for side in list('WESNR'):
        sect_list = seg[side]
        if len(sect_list) > 0:
            if side in ['W', 'S']:
                for sect in sect_list:
                    q_s, q_f, f_s, f_f, s_s, s_f = df.loc[sect,:]
                    if q_s > 0:
                        q[counter] = q_s
                        q_sal[counter] = 's'
                        Q[counter] = -q_f
                        Qf[counter] = -q_f
                        f[counter] = f_s
                        F[counter] = -f_f
                        s[counter] = s_s
                        S[counter] = s_f
                    elif q_s < 0:
                        q[counter] = q_f
                        Q[counter] = -q_s
                        f[counter] = f_f
                        F[counter] = -f_s
                        s[counter] = s_f
                        S[counter] = s_s
                    counter += 1
            elif side in ['E', 'N']:
                for sect in sect_list:
                    q_s, q_f, f_s, f_f, s_s, s_f = df.loc[sect,:]
                    if q_s > 0:
                        q[counter] =-q_f
                        Q[counter] = q_s
                        f[counter] = -f_f
                        F[counter] = f_s
                        s[counter] = s_f
                        S[counter] = s_s
                        
                    elif q_s < 0:
                        q[counter] = -q_s
                        q_sal[counter] = 's'
                        Q[counter] = q_f
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
                f[counter] = 0
                F[counter] = 0
    
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
        S = F/Q
    
        if testing == True:    
            print('\n' + seg_name + ':')
            print(' - Volume Flux adjustment = %0.5f' % (dq/Q.sum()))
            print(' -   Salt Flux adjustment = %0.5f' % (dqs/F.sum()))

    # The linear system of equations we will solve looks like
    # C * x = y => x = C-1 * y
    # where x is the vector of exchange fractions we are solving for,
    # C is the matrix of incoming fluxes,
    # and y is the vector of outgoing fluxes
    
    if N == 1: # H8 case, only a single section, no rivers
        A = np.array([1], ndmin=2)
        
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
            
            
        if (N == 3) and has_riv:
            # iterative method to scan over all possible river distributions
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
            if testing == True:
                print(good_guess)
            A = AA[good_guess,:,:].mean(axis=0)
        else:
            A = get_A(N, q, Q, f, F, yorig, A)
        
        if (N == 2) and (not has_riv):
            # check for bad solutions
            if (A < 0).any():
                print('Bad value for N=2 segment ' + seg_name)
                print('ORIGINAL')
                np.set_printoptions(precision=3, suppress=True)
                print(A)
                print('NEW')
                A = np.array([[0,1],[1,0]])
                np.set_printoptions(precision=3, suppress=True)
                print(A)
            
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

    if testing == True:
        print('Solution matrix A:')
        np.set_printoptions(precision=3, suppress=True)
        print(A)
        print(q_sal)
        
        print('')
        print(df_up_down)
    
    # cvec = 'rbgmc'
    # plt.close('all')
    # fig = plt.figure(figsize=(10,7))
    #
    #
    # plt.show()
    

