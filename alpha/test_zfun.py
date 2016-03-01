# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 08:45:10 2016

@author: PM5

Code to test the functions in the module zfun.

This is intended to be a chance to exercise the functions and look at them
easily, as opposed to a strict unit test.

"""
import numpy as np

from importlib import reload
import zfun
reload(zfun)

# choose the function to test
print('\n%s\n' % '** Choose function **')
fn_list_raw = dir(zfun)
fn_list = []
for fn in fn_list_raw:
    if fn[0] != '_':
        fn_list.append(fn)
Nfn = len(fn_list)
fn_dict = dict(zip(range(Nfn), fn_list))
for nfn in range(Nfn):
    print(str(nfn) + ': ' + fn_list[nfn])
my_nfn = int(input('-- Input number -- '))
fn_name = fn_dict[my_nfn]
fn = getattr(zfun, fn_name)

print('\n****** Info on function zfun.' + fn_name + ' ******\n')
print(help(fn))
print(50*'*' + '\n')

if fn_name == 'get_interpolant_fast':    
    xvec = np.linspace(0,10,11)
    # testing with some out of range values
    x = np.array([-1, .5, 12])   
    a = fn(x, xvec) 
    # and A is the correct answer
    A = np.array([[  0. ,   1. ,   0. ],
                  [  0. ,   1. ,   0.5],
                  [  9. ,  10. ,   1. ]])
    if (a==A).all():
        print('Success')
    else:
        print('Fail')
    

