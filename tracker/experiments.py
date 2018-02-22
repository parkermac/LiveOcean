#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 14:08:29 2018

@author: pm7
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
        
    elif ic_name == 'est1': # for the idealized estuary, all depths
        lonvec = np.linspace(-.8, 0, 20)
        latvec = np.linspace(44.8, 45.8, 40)
        lonmat, latmat = np.meshgrid(lonvec, latvec)
        plon_vec = lonmat.flatten()
        plat_vec = latmat.flatten()
        pcs_vec = np.arange(-.95, -.5, 10)
        
    return (gtagex, ic_name, plon_vec, plat_vec, pcs_vec)