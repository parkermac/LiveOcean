# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 16:01:15 2018

@author: lorenz
"""

import numpy as np

"""
input
x = Q(S)
comp = size of the window as an integer number
"""

def find_extrema(x,comp,min_transport):#new, reduced, better described
    #print(comp,min_transport)
    if min_transport<=10**(-10): #to kick out numerical errors!
        min_transport=10**(-10)
    #print(len(x),min_transport)
    #print(x[0])
    indices = []
    minmax = []
    i = 0
    while i < np.shape(x)[0]:
        if i-comp < 0:
            a = 0
        else:
            a=i-comp
        if i+comp>=len(x):
            b=None
        else:
            b=i+comp+1
        if x[i] == np.max(x[a:b]) and np.max(x[a:b]) != np.min(x[a:b]) and x[i] != x[i-1]:
            indices.append(i)
            minmax.append('max')
        elif x[i] == np.min(x[a:b]) and np.max(x[a:b]) != np.min(x[a:b]) and x[i] != x[i-1]:
            indices.append(i)
            minmax.append('min')
        i+=1
    #print(indices,minmax)
    #correct min min min or max max max parts, especially in the beginning and end of the Q(S)
    #index=[]
    ii=1
    while ii < len(indices):
        #print('minmin/maxmax',ii,len(indices))
        #print(minmax)
        index=[]
        if minmax[ii] == minmax[ii-1]:
            if minmax[ii] == 'max': #note the index of the smaller maximum
                if x[indices[ii]]>=x[indices[ii-1]]:
                    index.append(ii-1)
                else:
                    index.append(ii)
            elif minmax[ii] == 'min': #note the index of the greater minimum
                if x[indices[ii]]<=x[indices[ii-1]]:
                    index.append(ii-1)
                else:
                    index.append(ii)
            minmax = np.asarray(minmax)
            indices = np.asarray(indices)
            indices = np.delete(indices, index)
            minmax = np.delete(minmax, index)
        else:
            ii+=1
    #print(indices,minmax)
    #delete too small transports
    #print(indices,minmax)
    ii=0
    while ii < len(indices)-1: 
        index=[]
        #print('min_trans',ii,len(indices))
        #print(indices,minmax)
        if np.abs(np.abs(x[indices[ii+1]])-np.abs(x[indices[ii]])) <= min_transport:
            #print(np.abs(np.abs(x[indices[ii+1]])-np.abs(x[indices[ii]])),min_transport)
            if ii == 0: #if smin is involved and the transport is too small, smin has to change its min or max property
                #print('if')
                index.append(ii+1)
                if minmax[ii] == 'min':
                    minmax[ii] = 'max'
                else:
                    minmax[ii] = 'min'
            elif ii+1==len(indices)-1:#if smax is involved and the transport is too small, smin has to change its min or max property
                #print('elif')
                index.append(ii)
                if minmax[ii+1] == 'min':
                    minmax[ii+1] = 'max'
                else:
                    minmax[ii+1] = 'min'
            else: #else both involved div sals are kicked out
                #print('else')
                if ii+2 < len(indices)-1:
                #check and compare to i+2
                    if minmax[ii]=='min':
                        if x[indices[ii+2]]>x[indices[ii]]:
                            index.append(ii+2)
                            index.append(ii+1)
                        else:
                            index.append(ii)
                            index.append(ii+1)
                    elif minmax[ii]=='max':
                        if x[indices[ii+2]]<x[indices[ii]]:
                            index.append(ii+2)
                            index.append(ii+1)
                        else:
                            index.append(ii)
                            index.append(ii+1)
                else:
                    index.append(ii)
                    index.append(ii+1)
            #print(index)
            indices = np.delete(indices, index)
            minmax = np.delete(minmax, index)
            #print('after delete',indices,minmax)
        else:
            ii+=1
    #so far the first and last minmax does not correspond to smin and smax of the data, expecially smin due to numerical errors
    #correct smin index
    if len(x)>4:
        ii=1
        while np.abs(np.abs(x[ii])-np.abs(x[0])) < 10**(-10) and ii < len(x)-1:
            ii+=1
        indices[0]=ii-1
        #correct smax index
        if x[-1]==0: #for low salinity classes Q[-1] might not be zero as supposed.
            jj=-1
            while x[jj] == 0 and np.abs(jj) < len(x)-1:
                jj-=1
            indices[-1]=len(x)+jj+1
    return indices,minmax



"""
input
s=salinity array
Qv=Q(S)
Qs=Q^s(S)
comp=size of the comparison interval-window as an integer number
min_trans=minimum transport to consider
"""
def calc_bulk_values(s,Qv,Qs,comp,min_trans):
    smin=s[0]
    DeltaS=s[1]-s[0]
    ind,minmax = find_extrema(Qv,comp,min_trans) #use the find_extrema lgorithm
    #print(ind,minmax)
    div_sal=[]
    for i in range(0,len(ind)):#compute dividing salinities
            #print(Qvl[ind[i]])
        div_sal.append(smin+DeltaS*ind[i])
            #print(smin+dss*ind[i])
        #calculate transports etc.
    Q_in_m=[]
    Q_out_m=[]
    s_in_m=[]
    s_out_m=[]
    if minmax[0]=='max':#if inverse estuary, first outflow is set to zero
        Q_out_m.append(0)
        s_out_m.append(0)
    index=[]
    i=0
    while i < len(ind)-1:#compute the transports and sort to in and out
        Q_i=-(Qv[ind[i+1]]-Qv[ind[i]])
        F_i=-(Qs[ind[i+1]]-Qs[ind[i]])
        s_i=F_i/Q_i
        if Q_i<0 and np.abs(Q_i)>1:
            Q_out_m.append(Q_i)
            s_out_m.append(s_i)
        elif Q_i > 0 and np.abs(Q_i)>1:
            Q_in_m.append(Q_i)
            s_in_m.append(s_i)
        else:
            index.append(i)
        i+=1
    div_sal = np.delete(div_sal, index)
    # print(div_sal)
    # print(Q_in_m)
    # print(s_in_m)
    # print(Q_out_m)
    # print(s_out_m)
    return(Q_in_m,Q_out_m,s_in_m,s_out_m,div_sal)