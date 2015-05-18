"""
Plot fields in one or more history files.

On fjord this needs to be run with an X window.
"""

# setup
import os; import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(alp)
from datetime import datetime, timedelta
coast_file = Ldir['data'] + 'coast/pnw_coast_combined.mat'
import pfun; reload(pfun) # local module of plot routines

# choose the type of list to make
print 30*'*' + ' pan_plot ' + 30*'*'
print '\n%s\n' % '** Choose List type **'
lt_list = ['hand_selection','hindcast', 'mixed_hindcast', 'forecast', 'pnwtox']
Nlt = len(lt_list)
lt_dict = dict(zip(range(Nlt), lt_list))
for nlt in range(Nlt):
    print str(nlt) + ': ' + lt_list[nlt]
my_nlt = int(raw_input('-- Input number -- '))
list_type = lt_dict[my_nlt]

# choose the type of plot to make
print '\n%s\n' % '** Choose Plot type **'
pt_list_raw = dir(pfun)
pt_list = []
for pt in pt_list_raw:
    if pt[0] != '_':
        pt_list.append(pt)
Npt = len(pt_list)
pt_dict = dict(zip(range(Npt), pt_list))
for npt in range(Npt):
    print str(npt) + ': ' + pt_list[npt]
my_npt = int(raw_input('-- Input number -- '))
plot_type = pt_dict[my_npt]
whichplot = getattr(pfun, plot_type)

def make_fn_list(dt0,dt1,Ldir): # a helpful function
    from datetime import timedelta
    del_dt = timedelta(1)
    date_list = []
    fn_list = []
    dt = dt0
    while dt <= dt1:
        date_list.append(dt.strftime('%Y.%m.%d'))
        dt = dt + del_dt
    for dl in date_list:    
        f_string = 'f' + dl
        for nhis in range(2, 26): # range(2, 26) for a typical forecast              
            nhiss = ('0000' + str(nhis))[-4:]                        
            fn = (Ldir['roms'] + 'output/' + Ldir['gtag'] + '/'
                + f_string + '/ocean_his_' + nhiss + '.nc')
            fn_list.append(fn)
    return fn_list

# choose which file(s) to plot
fn_list = []
if list_type == 'hand_selection':
    # select one or more files using a dialog box
    from PySide import QtGui
    # unfortunately this kills the kernel in Anaconda
    fn_list, _ = QtGui.QFileDialog.getOpenFileNames(None, 'Choose File(s)',
        Ldir['roms'] + 'output/' + Ldir['gtag'] + '/')
elif list_type == 'hindcast':
    dt0 = datetime(2015,4,15) # first day
    dt1 = datetime(2015,4,15) # last day
    fn_list = make_fn_list(dt0,dt1,Ldir)
elif list_type == 'mixed_hindcast':
    fn_list1 = []
    fn_list2 = []
    dt0 = datetime(2013,1,1) # first day
    dt1 = datetime(2013,1,5) # last day
    fn_list = make_fn_list(dt0,dt1,Ldir)
    dt0 = datetime(2013,2,15) # first day
    dt1 = datetime(2013,2,19) # last day
    fn_list1 = make_fn_list(dt0,dt1,Ldir)
    for fn in fn_list1:
        fn_list.append(fn)
    dt0 = datetime(2013,3,15) # first day
    dt1 = datetime(2013,3,19) # last day
    fn_list2 = make_fn_list(dt0,dt1,Ldir)
    for fn in fn_list2:
        fn_list.append(fn)
elif list_type == 'forecast':
    ndays = 7
    dt0 = datetime.now() - timedelta(ndays) # first day
    dt1 = datetime.now() # last day
    fn_list = make_fn_list(dt0,dt1,Ldir)                    
elif list_type == 'pnwtox':
    # make a pnwtox-style list of files
    for nhis in range(1800, 1824): # have 1800 to 1929
        nhiss = ('0000' + str(nhis))[-4:]       
        fn = (Ldir['parent'] + 'roms/output/D2005_his/ocean_his_'
            + nhiss + '.nc')
        fn_list.append(fn)

# if plotting more than one file, prepare a directory for results
if len(fn_list) > 1:
    outdir0 = Ldir['LOo'] + 'plots/'
    Lfun.make_dir(outdir0, clean=False)
    outdir = outdir0 + list_type + '_' + plot_type + '/'
    Lfun.make_dir(outdir, clean=True)

# plot
if len(fn_list) == 1:
    # plot to screen
    whichplot(fn_list[0], alp, fn_coast=coast_file)
elif len(fn_list) > 1:
    # plot to a folder of files
    jj = 0     
    for fn in fn_list:    
        nouts = ('0000' + str(jj))[-4:]        
        outname = 'his_' + nouts + '.png'
        outfile = outdir + outname        
        print 'Plotting ' + fn       
        whichplot(fn, alp, fn_coast=coast_file,
            show_plot=False, save_plot=True, fn_out=outfile)
        jj += 1
        
