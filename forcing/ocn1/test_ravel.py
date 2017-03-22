#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 12:17:58 2017

@author: PM5

Test of raveling commands in numpy.

"""

import numpy as np

# set the size of the dimensions
L = 3
M = 3
N = 3

# create a data array (numbered layers)
a = np.array(range(N), dtype=float).reshape((N,1,1)) * np.ones((M,L))

# a flattened version of a
af = a.flatten()

# make an array the size of the array that contains a vertical index
b = np.ones((N,M,L), dtype=int)
b[1,1,1] = 2

bf = b.flatten()

LL = np.tile(np.arange(L), M*N)
MM = np.tile(np.repeat(np.arange(M), L), N)
NN = np.repeat(np.arange(N), M*L)

# this is the generic code for generating indices into the flattened
# array from indices in the original array
F = LL + MM*L + NN*M*L

# and this is indices into the flattened array at vertical levels given by b
f = LL + MM*L + bf*M*L

# get the flattened resuls
Af = af[f]

# And this is items in a at locations with vertical level b
A = Af.reshape(a.shape)


