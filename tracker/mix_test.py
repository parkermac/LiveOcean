"""
Code to test the random walk vertical mixing algorithm.
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand

h = 100
k0 = 0.01 # max value of k (m2 s-1)

def get_k(pz, k0, h):
    k = -4*k0*(pz/h)*((pz/h) + 1)
    return k
    
def get_dkdz(pz, k0, h):
    dkdz = -8*k0*pz/(h**2) - 4*k0/h
    return dkdz

# Result arrays packed as [time, particle]

NP = 1000 # number of particles
pz = np.linspace(-h+1,-1,NP)

NT = 10000 # number of timesteps
dt = 100 # timestep (s)
T = dt*np.arange(NT)
PT = T.reshape((NT,1)) * np.ones((1,NP))

PZ = np.nan * np.ones((NT, NP))
PZ[0,:] = pz

tt = 1
while tt < NT:

        pz0 = PZ[tt-1,:]
        k = get_k(pz0, k0, h)
        dkdz0 = get_dkdz(pz0, k0, h)
        pz_half = pz0 + 0.5*dt*dkdz0
        dkdz = get_dkdz(pz_half, k0, h)
        
        R = rand.randn(NP)
        pz1 = pz0 + dt*(R*np.sqrt(2*k/dt) + dkdz)
        
        # reflective boundaries
        pz1[pz1 < -h] = -2*h - pz1[pz1 < -h]
        pz1[pz1 > 0] = - pz1[pz1 > 0]
        
        
        PZ[tt,:] = pz1
        tt += 1

# plotting

plt.close('all')

fig = plt.figure(figsize=(14,8))

ax = fig.add_subplot(121)
ax.plot(PT[:,::50], PZ[:,::50], linewidth=.2)
ax.set_ylim(-h,0)
ax.set_xlim(0,NT*dt)

ax = fig.add_subplot(122)
ax.hist(pz0, bins=50, alpha=.3, orientation='horizontal')
ax.hist(pz1, bins=50, alpha=.3, orientation='horizontal')
ax.set_ylim(-h,0)

plt.show()