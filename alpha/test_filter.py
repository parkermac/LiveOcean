"""
Code to test the new zfun filtering functions for matrices.

Parker MacCready
"""

import numpy as np
import matplotlib.pyplot as plt

import zfun
from importlib import reload
reload(zfun)

# make some individual time series to filter
t = np.arange(100*24)
nt = len(t)
om = 2*np.pi/12
om1 = 2*np.pi/500
a1 = np.cos(om1*t) + .5*np.cos(om*t)
a2 = .7*np.sin(om1*t) + .5*np.cos(om*t)
a3 = .3*np.cos(om1*t) + .5*np.cos(om*t)
a4 = .5*np.sin(om1*t) + .5*np.cos(om*t)
# also pack them in an array with time along axis 0
A = np.nan * np.ones((nt,2,2))
A[:,0,0] = a1
A[:,1,0] = a2
A[:,0,1] = a3
A[:,1,1] = a4

to_test = 'godin'

if to_test == 'hanning':
    #  filter each one individually
    aa1 = zfun.filt_hanning(a1)
    aa2 = zfun.filt_hanning(a2)
    aa3 = zfun.filt_hanning(a3)
    aa4 = zfun.filt_hanning(a4)
    # and filter this using the function we are testing
    AA = zfun.filt_hanning_mat(A)
elif to_test == 'godin':
    #  filter each one individually
    aa1 = zfun.filt_godin(a1)
    aa2 = zfun.filt_godin(a2)
    aa3 = zfun.filt_godin(a3)
    aa4 = zfun.filt_godin(a4)
    # and filter this using the function we are testing
    AA = zfun.filt_godin_mat(A)

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(12,8))

ax = fig.add_subplot(111)
ax.plot(t,aa1,'-r',t,aa2, '-b',t,aa3, '-g',t,aa4, '-m', linewidth=10, alpha=.2)
ax.plot(t,AA[:,0,0],'-r',t,AA[:,1,0],'-b',t,AA[:,0,1],'-g',t,AA[:,1,1],'-m')
ax.set_title(to_test.title())
plt.show()