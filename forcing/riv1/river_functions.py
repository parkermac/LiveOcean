# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:11:16 2016

@author: PM5

Some extra river functions.

"""

def get_tc_rn(df):
    # makes a new column in df called 'tc_rn' which is the name of a
    # river that has T climatology

    for rn in df.index:
        if rn in ['coquille']:
            tc_rn = 'umpqua'
        elif rn in ['alsea']:
            tc_rn = 'siuslaw'
        elif rn in ['wilson', 'naselle', 'willapa', 'chehalis', 'humptulips',
                    'quinault', 'queets', 'hoh', 'calawah', 'hoko', 'elwha',
                    'dungeness', 'gold', 'sarita', 'sanjuan']:
            tc_rn = 'nehalem'
        elif rn in ['dosewallips', 'duckabush', 'hamma', 'skokomish', 'deschutes',
                    'nisqually', 'puyallup', 'green', 'snohomish',
                    'stillaguamish']:
            tc_rn = 'cedar'
        elif rn in ['skagit', 'samish']:
            tc_rn = 'nooksack'
        elif rn in ['clowhom', 'squamish']:
            tc_rn = 'fraser'
        elif rn in ['oyster', 'tsolum', 'englishman']:
            tc_rn = 'nanaimo'
        else:
            tc_rn = rn
        df.ix[rn, 'tc_rn'] = tc_rn

    return df
