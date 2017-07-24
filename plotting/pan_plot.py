"""
Plot fields in one or more history files.

On fjord this can be run from the command line, no X window needed,
but it is only for plotting to files, not the screen:

python pan_plot.py -x lobio3 -d 2013.01.02 -fno test.png -lt low_pass -pt P_basic

Running from the terminal on my mac, and making a movie:

python pan_plot.py -g aestus1 -t A1 -x ae1 -d 2013.02.07 -lt backfill -pt P_sectA -mov True
python pan_plot.py -g aestus1 -t A1 -x ae1 -d 2013.02.01 -lt backfill -pt P_sectA -mov True -nd 13

Running from the ipython command line:

cd /Users/PM5/Documents/LiveOcean/plotting
run pan_plot.py
run pan_plot.py -g aestus1 -t A1 -x ae1 -d 2013.02.07
run pan_plot.py -g cas1 -t f1 -x r820 -d 2013.01.01 -hs 25
run pan_plot.py -x lobio3 -d 2013.01.02 -fno test.png -lt low_pass -pt P_basic
run pan_plot.py -g aestus1 -t A1 -x ae1 -d 2013.02.07 -lt backfill -pt P_sectA -mov True
run pan_plot.py -g cascadia1 -t base -x lobio1 -d 2017.05.18 -lt snapshot -pt P_tracks
run pan_plot.py -d 2017.05.18 -lt merhab -pt P_tracks_MERHAB -mov True

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

#%% get optional command line arguments, any order
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', nargs='?', type=str,
                    default='cascadia1')
parser.add_argument('-t', '--tag', nargs='?', type=str,
                    default='base')
parser.add_argument('-x', '--ex_name', nargs='?', type=str,
                    default='lobio1')
parser.add_argument('-d', '--date_string', nargs='?', type=str,
                    default='2017.05.18')
parser.add_argument('-hs', '--hour_string', nargs='?', type=str,
                    default='02')
parser.add_argument('-nd', '--num_days', nargs='?', type=int,
                    default=0) # number of additional days
# more arguments that allow you to bypass the interactive choices
parser.add_argument('-lt', '--list_type', nargs='?', type=str,
                    default='')
parser.add_argument('-pt', '--plot_type', nargs='?', type=str,
                    default='')
parser.add_argument('-fno', '--fn_out', nargs='?', type=str,
                    default='')
parser.add_argument('-mov', '--make_movie', nargs='?', type=bool,
                    default=False)

args = parser.parse_args()

Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name

# choose the type of list to make
if len(args.list_type) == 0:
    print(30*'*' + ' pan_plot ' + 30*'*')
    print('\n%s\n' % '** Choose List type (return for snapshot) **')
    lt_list = ['snapshot', 'low_pass', 'backfill', 'forecast',
               'merhab', 'old_style',
               'atlantis', 'salish', 'salish_seq']
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
    fn_list.pop(0) # remove the first hour
elif list_type == 'backfill':
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
elif list_type=='salish':
    fn_list = []
    fn = (Ldir['parent'] + 'roms/output/salish_2006_4/ocean_his_5020.nc')
    fn_list.append(fn)
elif list_type=='salish_seq':
    fn_list = []
    for ii in range(4993, 5076): # have 4993 through 5075 on mac
        nstr = ('0000' + str(ii))[-4:]
        fn = (Ldir['parent'] + 'roms/output/salish_2006_4/ocean_his_' + nstr + '.nc')
        fn_list.append(fn)

#%% plot
in_dict = roms_plots.get_in_dict(plot_type)
vlims = in_dict['vlims']

if len(args.fn_out) == 0:
    if len(fn_list) == 1:
        # plot a single image to screen
        fn = fn_list[0]
        in_dict['fn'] = fn
        in_dict['fn_out'] = ''
        out_dict = whichplot(in_dict)
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
            in_dict['vlims'] = vlims
            out_dict = whichplot(in_dict)
            vlims = out_dict['vlims']
            jj += 1
        # and make a movie
        if args.make_movie:
            ff_str = ("ffmpeg -r 8 -pattern_type glob -i " + 
            outdir+"plot_%04d.png -vcodec libx264 -pix_fmt yuv420p -crf 25 "+outdir+"movie.mp4")
            os.system(ff_str)        
else:
    # plot a single image to a file
    fn = fn_list[0]
    in_dict['fn'] = fn
    in_dict['fn_out'] = args.fn_out
    out_dict = whichplot(in_dict)
 
