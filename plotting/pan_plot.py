"""
Plot fields in one or more history files.

On fjord this needs to be run with an X window.

Examples of running from the command line:

cd /Users/PM5/Documents/LiveOcean/plotting

run pan_plot.py -d 2015.02.01

run pan_plot.py -g cascadia2 -t frc2 -x lo1 -d 2013.01.09

run pan_plot.py -g cascadia2 -t frc2 -x lo1 -d 2013.01.01 -hs 25

"""

#%% setup
import os
import sys
import argparse
from datetime import datetime, timedelta
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
from importlib import reload
import Lfun
import roms_plots; reload(roms_plots)

#%% choices
in_dict = dict()

# COLOR LIMITS
vlims = dict()
# If you use () then the limits will be set by the first plot
# and then held constant at those levels thereafter.
vlims['salt'] = (28, 34)
vlims['temp'] = (8, 18)
vlims['NO3'] = (0, 40)
vlims['phytoplankton'] = (0,30)#(0, 40)
vlims['zooplankton'] = (0, 4)
vlims['oxygen'] = (0, 4) # for bottom DO (ml L-1)
vlims['TIC'] = (2000,2400)
vlims['alkalinity'] = (2000,2400)
in_dict['vlims'] = vlims

# OTHER
in_dict['z_level'] = -300 # z level to plot

#%% get optional command line arguments, any order
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', nargs='?', type=str,
                    default='cascadia1')
parser.add_argument('-t', '--tag', nargs='?', type=str,
                    default='base')
parser.add_argument('-x', '--ex_name', nargs='?', type=str,
                    default='lobio1')
parser.add_argument('-d', '--date_string', nargs='?', type=str,
                    default='2016.08.31')
parser.add_argument('-hs', '--hour_string', nargs='?', type=str,
                    default='02')
parser.add_argument('-nd', '--num_days', nargs='?', type=int,
                    default=0) # number of additional days
args = parser.parse_args()

Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name

# choose the type of list to make
print(30*'*' + ' pan_plot ' + 30*'*')
print('\n%s\n' % '** Choose List type (return for test) **')
lt_list = ['test', 'low_pass', 'hindcast', 'forecast', 'old_style', 'atlantis']
Nlt = len(lt_list)
lt_dict = dict(zip(range(Nlt), lt_list))
for nlt in range(Nlt):
    print(str(nlt) + ': ' + lt_list[nlt])

my_nlt = input('-- Input number -- ')
if len(my_nlt)==0:
    list_type = 'test'
else:
    list_type = lt_dict[int(my_nlt)]

dt0 = datetime.strptime(args.date_string, '%Y.%m.%d')
dt1 = dt0 + timedelta(args.num_days)

#%% choose the type of plot to make
print('\n%s\n' % '** Choose Plot type (return for P_basic) **')
pt_list_raw = dir(roms_plots)
pt_list = []
for pt in pt_list_raw:
    if pt[:2] == 'P_':
        pt_list.append(pt)
Npt = len(pt_list)
pt_dict = dict(zip(range(Npt), pt_list))
for npt in range(Npt):
    print(str(npt) + ': ' + pt_list[npt])
my_npt = input('-- Input number -- ')
if len(my_npt)==0:
    plot_type = 'P_basic'
else:
    plot_type = pt_dict[int(my_npt)]
whichplot = getattr(roms_plots, plot_type)

def make_fn_list(dt0, dt1, Ldir, hourmax=24):
    # a helpful function for making file lists
    from datetime import timedelta
    date_list = []
    fn_list = []
    dt = dt0
    while dt <= dt1:
        date_list.append(dt.strftime('%Y.%m.%d'))
        dt = dt + timedelta(1)
    for dl in date_list:
        f_string = 'f' + dl
        for nhis in range(2, hourmax+2):
            # range(2, 26) for a typical forecast
            nhiss = ('0000' + str(nhis))[-4:]
            fn = (Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' +
                  f_string + '/ocean_his_' + nhiss + '.nc')
            fn_list.append(fn)
    return fn_list

#%% choose which file(s) to plot
if list_type == 'test':
    # return a single default file name in the list
    fn_list = [Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' +
               'f' + args.date_string +
               '/ocean_his_00' + args.hour_string + '.nc']
elif list_type == 'hindcast':
    fn_list = make_fn_list(dt0,dt1,Ldir)
elif list_type == 'forecast':
    dt0 = datetime.now()
    fn_list = make_fn_list(dt0, dt0, Ldir, hourmax=72)
elif list_type == 'low_pass':
    date_list = []
    fn_list = []
    dt = dt0
    while dt <= dt1:
        date_list.append(dt.strftime('%Y.%m.%d'))
        dt = dt + timedelta(1)
    for dl in date_list:
        f_string = 'f' + dl
        fn = (Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
            + f_string + '/low_passed.nc')
        fn_list.append(fn)
elif list_type == 'old_style':
    # make a pnwtox-style list of files
    for nhis in range(1800, 1824): # have 1800 to 1929
        nhiss = ('0000' + str(nhis))[-4:]
        fn = (Ldir['parent'] + 'roms/output/D2005_his/ocean_his_' +
              nhiss + '.nc')
        fn_list.append(fn)
elif list_type=='atlantis':
    fn_list = []
    fn = (Ldir['parent'] + 'roms/output/salish_2006_4_lp/f2006.07.30/' +
        'low_passed.nc')
    fn_list.append(fn)

#%% plot
if len(fn_list) == 1:
    # plot to screen
    fn = fn_list[0]
    in_dict['fn'] = fn
    in_dict['fn_out'] = ''
    out_dict = whichplot(in_dict)
elif len(fn_list) > 1:
    #prepare a directory for results
    outdir0 = Ldir['LOo'] + 'plots/'
    Lfun.make_dir(outdir0, clean=False)
    outdir = outdir0 + list_type + '_' + plot_type + '_' + Ldir['gtagex'] + '/'
    Lfun.make_dir(outdir, clean=True)
    # plot to a folder of files
    jj = 0
    for fn in fn_list:
        nouts = ('0000' + str(jj))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = outdir + outname
        print('Plotting ' + fn)
        in_dict['fn'] = fn
        in_dict['fn_out'] = outfile
        in_dict['vlims'] = vlims
        out_dict = whichplot(in_dict)
        vlims = out_dict['vlims']
        jj += 1
