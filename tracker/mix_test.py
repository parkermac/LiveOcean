"""
Code to test the random walk vertical mixing algorithm.
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand

h = 20
k0 = 0.01 # max value of k (m2 s-1)

# warning: these functions don't check for out of bounds pz

def get_k(pz, k0, h):
    k = -4*k0*(pz/h)*((pz/h) + 1)
    return k
    
def get_dkdz(pz, k0, h):
    dkdz = -8*k0*pz/(h**2) - 4*k0/h
    return dkdz

# Result arrays packed as [time, particle]

NP = 1000 # number of particles
pz = np.linspace(-h+1,-1,NP)

Dmax = 100 # days to integrate
Tmax = 86400 * Dmax
dt = 3600 # timestep (s)
NT = int(Tmax/dt) # number of timesteps
T = dt*np.arange(NT)
PT = T.reshape((NT,1)) * np.ones((1,NP))

PZ = np.nan * np.ones((NT, NP))
PZ[0,:] = pz

tt = 1
while tt < NT:

        pz0 = PZ[tt-1,:]
        
        
        dkdz = get_dkdz(pz0, k0, h)
        
        pz_half = pz0 + 0.5*dt*dkdz
        
        k = get_k(pz_half, k0, h)
        
        R = rand.randn(NP)
        pz1 = pz0 + dt*(R*np.sqrt(2*k/dt) + dkdz)
        
        # reflective boundaries
        pz1[pz1 < -h] = -2*h - pz1[pz1 < -h]
        pz1[pz1 > 0] = - pz1[pz1 > 0]
        
        
        PZ[tt,:] = pz1
        tt += 1

# plotting

plt.close('all')

fig = plt.figure(figsize=(12,7))

ax = fig.add_subplot(121)
ssamp = 50
ax.plot(PT[:,::ssamp]/86400, PZ[:,::ssamp], linewidth=.5)
ax.set_ylim(-h,0)
dmax = NT*dt/86400
dmin = dmax - 20*dt/86400
ax.set_xlim(dmin, dmax)
ax.set_xlabel('Days')
ax.set_ylabel('Z (m)')
ax.set_title('Selected tracks, %d out of %d' % (int(NP/ssamp), NP))

ax = fig.add_subplot(122)
ax.hist(pz0, bins=50, alpha=.3, orientation='horizontal')
ax.hist(pz1, bins=50, alpha=.3, orientation='horizontal')
ax.set_ylim(-h,0)
ax.set_xlabel('Count')
ax.set_title('Histograms of last two timesteps, dt = %d s' % (dt))

plt.show()