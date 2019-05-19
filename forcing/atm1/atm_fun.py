"""
Functions to use with atmospheric forcing.  Translated from matlab to python.
"""

import numpy as np


invar_list = ['Q2', 'T2', 'PSFC', 'U10', 'V10','RAINCV', 'RAINNCV', 'SWDOWN', 'GLW']

outvar_list = ['Pair','rain','swrad','lwrad_down','Tair','Qair','Uwind','Vwind']

# for plotting, paired with outvar_list
lim_list = [(800,1050), (0,1e-6), (0,900), (0,400), (0,20), (50,100), (-10, 10), (-10,10)]

longname_dict = dict()
units_dict = dict()
timename_dict = dict()
for vn in outvar_list:
    if vn == 'Pair':
        nclongname = 'surface air pressure' 
        ncunits = 'millibar' 
        nctimename = 'pair_time' 
    elif vn == 'rain':
        nclongname = 'rain fall rate' 
        ncunits = 'kilograms meter-2 second-1' 
        nctimename = 'rain_time' 
    elif vn == 'swrad':
        nclongname = 'solar shortwave radiation flux' 
        ncunits = 'watts meter-2' 
        nctimename = 'srf_time' 
    elif vn == 'lwrad_down':
        nclongname = 'downwelling longwave radiation flux' 
        ncunits = 'watts meter-2' 
        nctimename = 'lrf_time' 
    elif vn == 'Tair':
        nclongname = 'surface air temperature' 
        ncunits = 'Celsius' 
        nctimename = 'tair_time' 
    elif vn == 'Qair':
        nclongname = 'surface air relative humidity' 
        ncunits = 'percentage' 
        nctimename = 'qair_time' 
    elif vn == 'Uwind':
        nclongname = 'surface u-wind component' 
        ncunits = 'meter second-1' 
        nctimename = 'wind_time' 
    elif vn == 'Vwind':
        nclongname = 'surface v-wind component' 
        ncunits = 'meter second-1' 
        nctimename = 'wind_time'
    longname_dict[vn] = nclongname
    units_dict[vn] = ncunits
    timename_dict[vn] = nctimename

def Z_wmo_RH(P,T,Q):
    # 5/21/2011 Nick Lederer, modified by Parker MacCready, and recoded
    # from matlab to python by PM 2019.05.16.  Tested against the matlab version
    # using Z_wmo_RH(1000, 10, .001) and both give 13.031628710406915.
    # 
    #  this converts mixing ratio (kg kg-1) which is the usual WRF output [CHECK!], into
    #  relative humidity (%) which is what ROMS expects
    # 
    #  INPUT:
    #  P in hectaPascal or millibar
    #  T in Celcius
    #  Q in kg kg-1
    # 
    #  OUTPUT:
    #  RH in percent
    #
    #  all equations come from Chapter 4 of
    #  http://www.wmo.int/pages/prog/www/IMOP/publications/CIMO-Guide/
    #  CIMO_Guide-7th_Edition-2008.html
    # 
    #  WMO GUIDE TO METEOROLOGICAL INSTRUMENTS AND METHODS OF OBSERVATION
    #           WMO-No. 8 (Seventh edition) (6 August 2008)
    #
    # Note 2019.05.15: this document no longer exists on the web.

    e_prime = Q*P/(0.62198+Q) # WHO equation 4.A.6
    fp = 1.0016 + 3.15e-6*P - 0.074/P # from Annex 4.B
    ew = 6.112*np.exp(17.62*T/(243.12 + T)) # from Annex 4.B
    ew_prime = fp*ew # from Annex 4.B
    RH = 100 * e_prime/ew_prime # from Annex 4.B

    return RH
    
    
if __name__ == '__main__':
    # examples of uses of the functions, executed in ipython as:
    # run atm_fun.py

    RH = Z_wmo_RH(1000, 10, .001)
    RH_expect = 13.031628710406915
    print('\nTest of Z_wmo_RH:')
    print(' Value - Expected Value = %0.1g\n' % (RH_expect - RH))
    
    print('\nTest of variable attributes')
    for vn in outvar_list:
        print('\n' + vn + ':')
        print(' longname = %s' % (longname_dict[vn]))
        print(' units = %s' % (units_dict[vn]))
        print(' timename = %s' % (timename_dict[vn]))
        