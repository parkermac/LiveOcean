"""
Plot fields in one or more history files.

Examples:

Plot a single figure to the screen with the default arguments
run pan_plot.py

Save multiple plots with color limits all set to match those set by
auto_lims() from the first plot:
run pan_plot.py -0 2017.09.01 -1 2017.09.03 -lt hourly -mov True

Use -avl False to have color limits all set to match those set by
pinfo.vlims_dict:

Example of an analytical run:
run pan_plot.py -g aestus1 -t A1 -x ae1 -0 2013.03.01

And a test of the MERHAB tracks:
run pan_plot.py -g cascadia1 -t base -x lobio5 -0 2013.01.31 -lt merhab -pt P_tracks_MERHAB -mov True

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

parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas4')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v2')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo6biom')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2018.09.29')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='')
# arguments that allow you to bypass the interactive choices
parser.add_argument('-hn', '--his_num', nargs='?', type=int, default=1)
parser.add_argument('-lt', '--list_type', nargs='?', type=str, default='')
parser.add_argument('-pt', '--plot_type', nargs='?', type=str, default='')
# arguments that influence other behavior
#  e.g. make a movie, override auto color limits
parser.add_argument('-mov', '--make_movie', default=False, type=boolean_string)
parser.add_argument('-avl', '--auto_vlims', default=True, type=boolean_string)
parser.add_argument('-test', '--testing', default=False, type=boolean_string)

args = parser.parse_args()
if len(args.date_string1) == 0:
    args.date_string1 = args.date_string0

Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name

# choose the type of list to make
if len(args.list_type) == 0:
    print(30*'*' + ' pan_plot ' + 30*'*')
    print('\n%s\n' % '** Choose List type (return for snapshot) **')
    lt_list = ['snapshot', 'low_passed', 'daily', 'hourly ', 'forecast',
               'merhab', 'allhours']
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

# choose the type of plot to make
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
if plot_type == 'P_tracks_MERHAB':
    # enforce list_type
    if list_type != 'merhab':
        list_type = 'merhab'
        print('NOTE: Overriding chosen list_type and using merhab instead.')

# get list of history files to plot
fn_list = Lfun.get_fn_list(list_type, Ldir,
    args.date_string0, args.date_string1, his_num=args.his_num)

if (list_type == 'merhab') and (args.testing == True):
    fn_list = fn_list[:10]
    
if (list_type == 'forecast') and (args.testing == True):
    # list of all history files in a directory
    hourmax = 4 #Ldir['forecast_days'] * 24
    dt0 = datetime.strptime(args.date_string0, '%Y.%m.%d')
    fn_list = Lfun.fn_list_utility(dt0, dt0, Ldir, hourmax=hourmax)
    
# PLOTTING
in_dict = dict()
in_dict['auto_vlims'] = args.auto_vlims
in_dict['testing'] = args.testing

if plot_type == 'P_3day':
    # special additions for this plot type
    in_dict['list_type'] = list_type
    fn_list = [fn_list[0], fn_list[0]] # two repeats so that we save to a file

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
