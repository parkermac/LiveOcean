#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 14:08:29 2018

@author: pm7

This is where you set the run "gtagex" and the initial condition
based on an experiment name passed by the calling code.

"""

import numpy as np

def make_ic(exp_name):
    
    if exp_name == 'jdf5':
        gtagex = 'cascadia1_base_lobio5'
        ic_name = 'jdf0'
    elif exp_name == 'jdf6':
        gtagex = 'cas3_v0_lo6m'
        ic_name = 'jdf0'
    elif exp_name == 'hc5':
        gtagex = 'cascadia1_base_lobio5'
        ic_name = 'hc0'
    elif exp_name == 'hc6':
        gtagex = 'cas3_v0_lo6m'
        ic_name = 'hc0'
    elif exp_name == 'ae0':
        gtagex = 'aestus1_A1_ae1'
        ic_name = 'est0'
    elif exp_name == 'ae1':
        gtagex = 'aestus1_A1_ae1'
        ic_name = 'est1'
    elif exp_name == 'ae2':
        gtagex = 'aestus1_A1_ae1'
        ic_name = 'est2'
    elif exp_name == 'ae3':
        gtagex = 'aestus1_A1_ae1'
        ic_name = 'est3'
        
    # routines to set particle initial locations, all numpy arrays
    #
    # first create three vectors of initial locations
    # plat00 and plon00 should be the same length,
    # and the length of pcs00 is however many vertical positions you have at
    # each lat, lon (expressed as fraction of depth -1 < pcs < 1)
    #
    if ic_name == 'hc0': # Hood Canal
        lonvec = np.linspace(-122.65, -122.45, 30)
        latvec = np.linspace(47.2, 47.35, 30)
        lonmat, latmat = np.meshgrid(lonvec, latvec)
        plon_vec = lonmat.flatten()
        plat_vec = latmat.flatten()
        pcs_vec = np.array([-.05])
        
    elif ic_name == 'jdf0': # Mid-Juan de Fuca
        lonvec = np.linspace(-123.85, -123.6, 20)
        latvec = np.linspace(48.2, 48.4, 20)
        lonmat, latmat = np.meshgrid(lonvec, latvec)
        plon_vec = lonmat.flatten()
        plat_vec = latmat.flatten()
        pcs_vec = np.array([-.05])
        
    elif ic_name == 'est0': # for the idealized estuary, surface
        lonvec = np.linspace(-.5, .5, 20)
        latvec = np.linspace(44.7, 45.3, 20)
        lonmat, latmat = np.meshgrid(lonvec, latvec)
        plon_vec = lonmat.flatten()
        plat_vec = latmat.flatten()
        pcs_vec = np.array([-.05])
        
    elif ic_name == 'est1': # for the idealized estuary
        lonvec = np.linspace(-.8, 0, 20)
        latvec = np.linspace(44.8, 45.8, 40)
        lonmat, latmat = np.meshgrid(lonvec, latvec)
        plon_vec = lonmat.flatten()
        plat_vec = latmat.flatten()
        pcs_vec = np.linspace(-.95, -.5, 10)
        
    elif ic_name == 'est2': # for the idealized estuary
        lonvec = np.linspace(-0.5, 0.5, 20)
        latvec = np.linspace(44.8, 45.2, 10)
        lonmat, latmat = np.meshgrid(lonvec, latvec)
        plon_vec = lonmat.flatten()
        plat_vec = latmat.flatten()
        pcs_vec = np.linspace(-.95, -.5, 10)
        
    elif ic_name == 'est3': # for the idealized estuary
        lonvec = np.linspace(-.2, .5, 5)
        latvec = np.linspace(44.9, 45.1, 10)
        lonmat, latmat = np.meshgrid(lonvec, latvec)
        plon_vec = lonmat.flatten()
        plat_vec = latmat.flatten()
        pcs_vec = np.linspace(-.95, -.1, 5)
    
    # Create full output vectors (each has one value per point).  This
    # code takes each lat, lon location and then assigns it to NSP points
    # corresponding to the vector of pcs values.  However you could write a
    # different version that only released points below a certain depth,
    # or other criterion.
    if len(plon_vec) != len(plat_vec):
        print('WARNING: Problem with length of initial lat, lon vectors')
    NSP = len(pcs_vec)
    NXYP = len(plon_vec)
    plon_arr = plon_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    plat_arr = plat_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    pcs_arr = np.ones((NXYP,NSP)) * pcs_vec.reshape(1,NSP)
    plon00 = plon_arr.flatten()
    plat00 = plat_arr.flatten()
    pcs00 = pcs_arr.flatten()
        
    return (gtagex, ic_name, plon00, plat00, pcs00)