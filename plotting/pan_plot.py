"""
Plot fields in one or more history files.

On fjord this needs to be run with an X window.

Rewritten 2012.12.02 to put more choices about which files to plot
on the command line, thereby avioding the use of tkinter, which crashed often.

"""

#%% setup
import os
import sys
import argparse
from datetime import datetime, timedelta
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import pfun
from importlib import reload
reload(pfun)

# get optional command line arguments, any order
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cascadia1')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='base')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo1')
parser.add_argument('-d', '--date_string', nargs='?', type=str, default='2015.09.19')
parser.add_argument('-hs', '--hour_string', nargs='?', type=str, default='02')
# num_days = number of additional days
parser.add_argument('-nd', '--num_days', nargs='?', type=int, default=0)
args = parser.parse_args()

Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name
coast_file = Ldir['data'] + 'coast/pnw_coast_combined.mat'

# choose the type of list to make
print(30*'*' + ' pan_plot ' + 30*'*')
print('\n%s\n' % '** Choose List type **')
lt_list = ['test', 'low_pass', 'hindcast', 'forecast', 'old_style']
Nlt = len(lt_list)
lt_dict = dict(zip(range(Nlt), lt_list))
for nlt in range(Nlt):
    print(str(nlt) + ': ' + lt_list[nlt])
my_nlt = int(input('-- Input number -- '))
list_type = lt_dict[my_nlt]

dt0 = datetime.strptime(args.date_string, '%Y.%m.%d')
dt1 = dt0 + timedelta(args.num_days)

#%% choose the type of plot to make
print('\n%s\n' % '** Choose Plot type **')
pt_list_raw = dir(pfun)
pt_list = []
for pt in pt_list_raw:
    if pt[0] != '_':
        pt_list.append(pt)
Npt = len(pt_list)
pt_dict = dict(zip(range(Npt), pt_list))
for npt in range(Npt):
    print(str(npt) + ': ' + pt_list[npt])
my_npt = int(input('-- Input number -- '))
plot_type = pt_dict[my_npt]
whichplot = getattr(pfun, plot_type)

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
        for nhis in range(2, hourmax+2):  # range(2, 26) for a typical forecast
            nhiss = ('0000' + str(nhis))[-4:]
            fn = (Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' +
                  f_string + '/ocean_his_' + nhiss + '.nc')
            fn_list.append(fn)
    return fn_list

#%% choose which file(s) to plot (way too complicated)
if list_type == 'test':
    # return a single default file name in the list
    fn_list = [Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' +
               'f' + args.date_string +
               '/ocean_his_00' + args.hour_string + '.nc']
#    fn_list = ['/Users/PM5/Documents/roms/output/salish_2006_4/ocean_his_8701.nc']
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

#%% if plotting more than one file, prepare a directory for results
if len(fn_list) > 1:
    outdir0 = Ldir['LOo'] + 'plots/'
    Lfun.make_dir(outdir0, clean=False)
    outdir = outdir0 + list_type + '_' + plot_type + '/'
    Lfun.make_dir(outdir, clean=True)

#%% plot
if len(fn_list) == 1:
    # plot to screen
    whichplot(fn_list[0], alp, Ldir, fn_coast=coast_file)
elif len(fn_list) > 1:
    # plot to a folder of files
    jj = 0
    for fn in fn_list:
        nouts = ('0000' + str(jj))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = outdir + outname
        print('Plotting ' + fn)
        whichplot(fn, alp, Ldir, fn_coast=coast_file,
            show_plot=False, save_plot=True, fn_out=outfile)
        jj += 1
