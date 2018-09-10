"""
Plots mooring records.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime

import pandas as pd

import os
import sys

pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import zfun

pth = os.path.abspath('../plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'extract/'
# choose the mooring extraction to plot

plt.close('all')
fs1=16 # fontsize for labels
fs2 = 12
for ww in ['0','1','2','3']:
    
    moor_file = ('moor_cas4_v1_lo6biom_willapa' + ww
        + '_hourly_2017.01.01_2017.12.31.nc')
    fn = indir + moor_file
    ds = nc.Dataset(fn)

    #%% load and organize data

    # create a dict into which to load everything
    V = dict()
    # and a dict of units (not used, but for reference)
    Vu = dict()

    # gather some fields
    vn_list = ['salt', 'temp',
        'NO3', 'phytoplankton',
        'TIC', 'alkalinity',
        'PH', 'ARAG', 'z_rho', 'svstr', 'ocean_time']

    # gather data (including ocean_time and units
    for vn in vn_list:
        V[vn] = ds[vn][:]
        try:
            Vu[vn] = ds[vn].units
        except AttributeError:
            Vu[vn] = ''

    lon = ds['lon_rho'][:]
    lat = ds['lat_rho'][:]
    
    NMID = 14
    ztop = V['z_rho'][-1,:].mean()
    zmid = V['z_rho'][NMID,:].mean()
    zbot = V['z_rho'][0,:].mean()

    ds.close()

    # load river files
    riv_list = ['columbia', 'willapa', 'naselle']
    riv_dict = dict()
    for riv in riv_list:
        riv_dict[riv] = pd.read_pickle(Ldir['parent'] +
            'ptools_output/river/coastal_2017/' + riv + '.p')

    #%% plotting
    fig = plt.figure(figsize=(20,8))

    mdays = Lfun.modtime_to_mdate_vec(V['ocean_time'])
    mdt = mdates.num2date(mdays) # list of datetimes of data

    nn_dict = {'salt': 1,
            'NO3': 2,
            'PH': 3,
            'temp': 5,
            'phytoplankton': 6,
            'ARAG': 7,
            'TIC': 4,
            'alkalinity': 8}
        
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

    for vn in nn_dict:
        ax = fig.add_subplot(3,4,nn_dict[vn])
        ax.plot(mdt, fac_dict[vn] * zfun.filt_godin(V[vn][0,:]), '-b')
        ax.plot(mdt, fac_dict[vn] * zfun.filt_godin(V[vn][NMID,:]), '-g')
        ax.plot(mdt, fac_dict[vn] * zfun.filt_godin(V[vn][-1,:]), '-r')
        ax.grid(True)
        ax.set_xlim(mdt[0], mdt[-1])
        ax.ticklabel_format(useOffset=False, axis='y')
        ax.text(.05, .85, tstr_dict[vn],
            horizontalalignment='left', transform=ax.transAxes, fontsize=fs1)
        ax.text(.05, .7, units_dict[vn],
            horizontalalignment='left', transform=ax.transAxes, fontsize=fs1)
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.set_xticklabels([])
        ax.set_ylim(vlims_dict[vn])
        # add depth labels
        if vn=='salt':
            ax.text(.95, .05, ('Z = %d (m)' % (int(zbot))),
                horizontalalignment='right', color='b',
                transform=ax.transAxes, fontsize=fs2, fontweight='bold')
            ax.text(.95, .15, ('Z = %d (m)' % (int(zmid))),
                horizontalalignment='right', color='g',
                transform=ax.transAxes, fontsize=fs2, fontweight='bold')
            ax.text(.95, .25, ('Z = %d (m)' % (int(ztop))),
                horizontalalignment='right', color='r',
                transform=ax.transAxes, fontsize=fs2, fontweight='bold')

    vn = 'svstr'
    ax = fig.add_subplot(3,4,9)
    fld = zfun.filt_AB8d(V[vn][:])
    ax.plot(mdt, fld, '-k')
    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages
    ax.fill_between(mdt, fld, where=(fld<=0), color='slateblue', alpha=0.5)
    ax.fill_between(mdt, fld, where=(fld>=0), color='lightsalmon', alpha=0.5)
    ax.grid(True)
    ax.set_xlim(mdt[0], mdt[-1])
    ax.ticklabel_format(useOffset=False, axis='y')
    ax.text(.05, .85, 'N-S Windstress (8 day filter)',
        horizontalalignment='left', transform=ax.transAxes, fontsize=fs1)
    ax.text(.05, .7, '$(Pa)$',
        horizontalalignment='left', transform=ax.transAxes, fontsize=fs1)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.set_xticklabels([])
    ax.set_ylim(-.1, .3)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.set_xlabel('Date 2017')
    ax.xaxis.set_tick_params(labelrotation=45)

    # rivers
    ax = fig.add_subplot(3,4,10)
    rr = riv_dict['columbia']
    rr = rr/1000
    rr.plot(style='-', color='brown')
    ax.grid(True)
    ax.set_xlim(mdt[0], mdt[-1])
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.set_ylim(0, 20)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.set_xlabel('Date 2017')
    ax.xaxis.set_tick_params(labelrotation=45)
    ax.text(.05, .85, 'Columbia River Flow $(1000\ m^{3}s^{-1})$',
        horizontalalignment='left', transform=ax.transAxes,
        color='brown', fontsize=fs1)

    # rivers
    ax = fig.add_subplot(3,4,11)
    rr = riv_dict['willapa']
    rr = rr/1000
    rr.plot(style='-', color='orange')
    rr = riv_dict['naselle']
    rr = rr/1000
    rr.plot(style='-', color='teal')
    ax.text(.05, .85, 'Willapa River Flow $(1000\ m^{3}s^{-1})$',
        horizontalalignment='left', transform=ax.transAxes,
        color='orange', fontsize=fs1)
    ax.text(.05, .7, 'Naselle River Flow $(1000\ m^{3}s^{-1})$',
        horizontalalignment='left', transform=ax.transAxes,
        color='teal', fontsize=fs1)
    ax.grid(True)
    ax.set_xlim(mdt[0], mdt[-1])
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.set_ylim(0, .5)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.set_xlabel('Date 2017')
    ax.xaxis.set_tick_params(labelrotation=45)

    # map
    ax = fig.add_subplot(3,4,12)
    pfun.add_coast(ax)
    ax.axis([-124.6, -123.6, 46.3, 46.8])
    pfun.dar(ax)
    ax.plot(lon, lat, '*r')
    sta = moor_file.split('_')[4]
    ax.set_title(sta)
    
    fig.tight_layout()
    
    out_dir = Ldir['LOo']+'willapa_moorings/'
    Lfun.make_dir(out_dir)
    plt.savefig(out_dir + sta + '.png')
    
    plt.close()

