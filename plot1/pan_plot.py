"""
Plot fields in one or more history files.

Examples:

Plot a single figure to the screen:
run pan_plot.py

Save multiple plots with color limits all set to match those set by
auto_lims() from the first plot:
run pan_plot.py -g cascadia1 -t base -x lobio5 -d 2013.01.31 -lt backfill

Save multiple plots with color limits all set to match those set by
pinfo.vlims_dict:
run pan_plot.py -g cascadia1 -t base -x lobio5 -d 2013.01.31 -lt backfill -avl False

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

def boolean_string(s):
    if s not in ['False', 'True']:
        raise ValueError('Not a valid boolean string')
    return s == 'True' # note use of ==

#%% get optional command line arguments, any order
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', nargs='?', type=str,
                    default='cas4')
parser.add_argument('-t', '--tag', nargs='?', type=str,
                    default='v1')
parser.add_argument('-x', '--ex_name', nargs='?', type=str,
                    default='lo6biom')
parser.add_argument('-d', '--date_string', nargs='?', type=str,
                    default='2017.05.08')
parser.add_argument('-hs', '--hour_string', nargs='?', type=str,
                    default='01')
parser.add_argument('-nd', '--num_days', nargs='?', type=int,
                    default=0) # number of additional days
# more arguments that allow you to bypass the interactive choices
parser.add_argument('-lt', '--list_type', nargs='?', type=str,
                    default='')
parser.add_argument('-pt', '--plot_type', nargs='?', type=str,
                    default='')
parser.add_argument('-mov', '--make_movie', default=False, type=boolean_string)
parser.add_argument('-avl', '--auto_vlims', default=True, type=boolean_string)

args = parser.parse_args()

Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name

# choose the type of list to make
if len(args.list_type) == 0:
    print(30*'*' + ' pan_plot ' + 30*'*')
    print('\n%s\n' % '** Choose List type (return for snapshot) **')
    lt_list = ['snapshot', 'low_pass', 'backfill', 'forecast',
               'merhab']
    Nlt = len(lt_list)
    lt_dict = dict(zip(range(Nlt), lt_list))
    for nlt in range(Nlt):
        print(str(nlt) + ': ' + lt_list[nlt])
    my_nlt = input('-- Input number -- ')
    if len(my_nlt)==0:
        list_type = 'snapshot'
    else:
        list_type = lt_dict[int(my_nlt)]
else:
    list_type = args.list_type

dt0 = datetime.strptime(args.date_string, '%Y.%m.%d')
dt1 = dt0 + timedelta(args.num_days)

#%% choose the type of plot to make
if len(args.plot_type) == 0:
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
else:
    plot_type = args.plot_type
    
whichplot = getattr(roms_plots, plot_type)

def make_fn_list(dt0, dt1, Ldir, hourmin=0, hourmax=24):
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
        for nhis in range(hourmin+1, hourmax+2):
            nhiss = ('0000' + str(nhis))[-4:]
            fn = (Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' +
                  f_string + '/ocean_his_' + nhiss + '.nc')
            fn_list.append(fn)
    return fn_list

#%% choose which file(s) to plot
if list_type == 'snapshot':
    # return a single default file name in the list
    fn_list = [Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' +
               'f' + args.date_string +
               '/ocean_his_00' + args.hour_string + '.nc']
elif plot_type == 'P_tracks_MERHAB':
    # enforce list_type
    if list_type != 'merhab':
        print('Need to use list_type=merhab for plot_type=P_tracks_MERHAB')
    # get a list of all but the first input file
    in_dir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' +
           'f' + args.date_string + '/')
    fn_list_raw = os.listdir(in_dir)
    fn_list = [(in_dir + ff) for ff in fn_list_raw if 'ocean_his' in ff]
    fn_list.sort()
    # testing - make a shorter list
    #fn_list = fn_list[:10]
    fn_list.pop(0) # remove the first hour
elif list_type == 'backfill':
    fn_list = make_fn_list(dt0,dt1,Ldir)
    fn_list = fn_list[:4] # testing
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
        
# plot

in_dict = dict()
in_dict['auto_vlims'] = args.auto_vlims

if len(fn_list) == 1:
    # plot a single image to screen
    fn = fn_list[0]
    in_dict['fn'] = fn
    in_dict['fn_out'] = ''
    whichplot(in_dict)
elif len(fn_list) > 1:
    # prepare a directory for results
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
        whichplot(in_dict)
        # after the first plot we no longer change vlims
        in_dict['auto_vlims'] = False
        jj += 1
    # and make a movie
    if args.make_movie:
        ff_str = ("ffmpeg -r 8 -i " + 
        outdir+"plot_%04d.png -vcodec libx264 -pix_fmt yuv420p -crf 25 "
        +outdir+"movie.mp4")
        os.system(ff_str)
