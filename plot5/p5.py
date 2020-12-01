"""
Plot fields in one or more history files.

"""

#%% setup
import os, sys
import argparse
from datetime import datetime, timedelta
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
import zrfun
import netCDF4 as nc
import numpy as np

from importlib import reload
import p5_plots
reload(p5_plots)

import matplotlib.pyplot as plt
plt.close('all')

parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', type=str, default='cas6')
parser.add_argument('-t', '--tag', type=str, default='v3')
parser.add_argument('-x', '--ex_name', type=str, default='lo8b')
parser.add_argument('-0', '--date_string0', type=str, default='2019.07.04')
parser.add_argument('-1', '--date_string1', type=str, default='')
# arguments that allow you to bypass the interactive choices
parser.add_argument('-hn', '--his_num', type=int, default=1)
parser.add_argument('-lt', '--list_type', type=str, default='')
parser.add_argument('-pt', '--plot_type', type=str, default='')
# arguments that influence other behavior
#  e.g. make a movie, override auto color limits
parser.add_argument('-mov', '--make_movie', default=False, type=zfun.boolean_string)
parser.add_argument('-avl', '--auto_vlims', default=True, type=zfun.boolean_string)
parser.add_argument('-test', '--testing', default=False, type=zfun.boolean_string)

args = parser.parse_args()
in_dict = args.__dict__

if len(in_dict['date_string1']) == 0:
    in_dict['date_string1'] = in_dict['date_string0']
    
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + in_dict['ex_name']

# choose the type of list to make
if len(args.list_type) == 0:
    print(30*'*' + ' pan_plot ' + 30*'*')
    print('\n%s\n' % '** Choose List type (return for snapshot) **')
    lt_list = ['snapshot', 'daily', 'hourly ', 'allhours']
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
    list_type = in_dict['list_type']

# choose the type of plot to make
if len(args.plot_type) == 0:
    print('\n%s\n' % '** Choose Plot type (return for P_salt_1) **')
    pt_list_raw = dir(p5_plots)
    pt_list = [item for item in pt_list_raw if item[:2] == 'P_']
    Npt = len(pt_list)
    pt_dict = dict(zip(range(Npt), pt_list))
    for npt in range(Npt):
        print(str(npt) + ': ' + pt_list[npt])
    my_npt = input('-- Input number -- ')
    if len(my_npt)==0:
        plot_type = 'P_salt_1'
    else:
        plot_type = pt_dict[int(my_npt)]
else:
    plot_type = in_dict['plot_type']
whichplot = getattr(p5_plots, plot_type)

# get list of history files to plot
fn_list = Lfun.get_fn_list(list_type, Ldir,
    in_dict['date_string0'], in_dict['date_string1'], his_num=in_dict['his_num'])
    
# get a mooring record (eventually put in a function)
m_fn_list = Lfun.get_fn_list(in_dict['list_type'], Ldir, in_dict['date_string0'], in_dict['date_string1'])
ot_list = []
zeta_list = []
uwind_list = []
vwind_list = []
G = zrfun.get_basic_info(m_fn_list[0], only_G=True)
m_lon = -124.5; m_lat = 47
mi = zfun.find_nearest_ind(G['lon_rho'][0,:], m_lon)
mj = zfun.find_nearest_ind(G['lat_rho'][:,0], m_lat)
for fn in m_fn_list:
    ds = nc.Dataset(fn)
    ot_list.append(ds['ocean_time'][0])
    zeta_list.append(ds['zeta'][0,mj,mi])
    uwind_list.append(ds['Uwind'][0,mj,mi])
    vwind_list.append(ds['Vwind'][0,mj,mi])
    ds.close()
ot_vec = zfun.fillit(np.array(ot_list))
zeta_vec = zfun.fillit(np.array(zeta_list))
uwind_vec = zfun.fillit(np.array(uwind_list))
vwind_vec = zfun.fillit(np.array(vwind_list))
in_dict['ot_vec'] = ot_vec
in_dict['zeta_vec'] = zeta_vec
in_dict['uwind_vec'] = uwind_vec
in_dict['vwind_vec'] = vwind_vec
in_dict['m_lon'] = m_lon
in_dict['m_lat'] = m_lat
        
# PLOTTING

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
