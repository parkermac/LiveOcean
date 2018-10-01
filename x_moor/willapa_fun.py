"""
Functions for plotting Willapa Moorings.
"""

vlims_dict = {'salt': (14, 35),
        'temp': (0, 24),
        'NO3': (0, 44),
        'phytoplankton': (0,45),
        'zooplankton': (0, 4),
        'oxygen': (0, 8),
        'TIC': (1400, 2400),
        'alkalinity': (1400,2400),
        'PH': (7, 8.5),
        'ARAG': (0, 3),
        'Ldetritus': ()}
        
# Units (after multiplying by scaling factor)
units_dict = {'salt': '',
             'temp': ' $(^{\circ}C)$',
             'NO3': ' $(\mu mol\ L^{-1})$',
             'phytoplankton': ' $(\mu g\ chl\ C\ L^{-1})$',
             'zooplankton': ' $(\mu g\ chl\ C\ L^{-1})$',
             'oxygen': ' $(ml\ L^{-1})$',
             'TIC': ' $(\mu mol\ kg^{-1})$',
             'alkalinity': ' $(\mu\ equivalents\ kg^{-1})$',
             'PH': '',
             'ARAG': ''}
             
# Scaling factors
fac_dict =  {'salt': 1,
             'temp': 1,
             'NO3': 1,
             'phytoplankton': 2.5,
             'zooplankton': 2.5,
             'oxygen': 0.032/1.42903, # convert mmol m-3 to ml L-1
             'TIC': 1000/1025, # convert L-1 to kg-1
             'alkalinity': 1000/1025,
             'PH': 1,
             'ARAG': 1}
             
# String form to use in titles
tstr_dict = {'salt': 'Salinity',
             'temp': 'Temperature',
             'NO3': 'Nitrate',
             'phytoplankton': 'Phytoplankton',
             'zooplankton': 'Zooplankton',
             'oxygen': 'DO',
             'TIC': 'DIC',
             'alkalinity': 'Alkalinity',
             'PH': 'pH',
             'ARAG': '$\Omega_{arag}$'}