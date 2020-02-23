"""
Plot the results of the flux age engine as a movie.

This is aimed to be used especially with the "S_" integrations,
meaning like a steady river or ocean source

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
import zfun

import tef_fun
import flux_fun

# select the indir
indir0 = Ldir['LOo'] + 'tef/flux_engine/cas6_v3_lo8b/'
indir = indir0 + 'flux/'

# load the DataFrame of results of flux_engine.py
infile = Lfun.choose_item(indir, tag='S_', exclude_tag='AGE',
    itext='Choose flux engine output file:')
aa = pd.read_pickle(indir + infile)

# create the outdir
outname = infile.replace('aa_','').replace('.p','')
outdir = indir + 'movie_' + outname +'/'
Lfun.make_dir(outdir, clean=True)


# load a Series of the volumes of each segment, created by flux_get_vol.py
v_df = pd.read_pickle(indir + 'volumes.p')

plt.close('all')

fs = 16
lw = 3

Dist_dict = {}
for ch in flux_fun.channel_dict.keys():
    seg_list = flux_fun.seg_dict[ch]
    dist = flux_fun.make_dist(v_df.loc[seg_list,'lon'],v_df.loc[seg_list,'lat'])
    Dist = pd.Series(dist, index=seg_list)
    Dist_dict[ch] = Dist
# hand edits
Dist_dict['Admiralty Inlet to South Sound'] += Dist_dict['Juan de Fuca to Strait of Georgia']['J4']
Dist_dict['Hood Canal'] += Dist_dict['Admiralty Inlet to South Sound']['A3']
Dist_dict['Whidbey Basin'] += Dist_dict['Admiralty Inlet to South Sound']['M1']

testing = True

tind = aa.index

if testing:
    tind = tind[:5]
else:
    pass
    #tind = tind[:51]
    

ii = 0
for tt in tind:
    
    fig = plt.figure(figsize=(13,8))
    ax1 = fig.add_subplot(111)

    ch_counter = 0
    for ch in flux_fun.channel_dict.keys():
        dist = Dist_dict[ch]
        seg_list = flux_fun.seg_dict[ch]

        vs = [s + '_s' for s in seg_list]
        vf = [s + '_f' for s in seg_list]

        color = flux_fun.clist[ch_counter]
        ax1.plot(dist, aa.loc[tt,vs].values,'-',color=color, linewidth=lw)
        ax1.plot(dist, aa.loc[tt,vf].values,'--',color=color, linewidth=lw)
        ch_counter += 1
    
    ax1.set_ylim(0,1.1)
    ax1.set_xlim(-10,410)

    ax1.grid(True)

    ax1.set_ylabel('Concentration', fontsize=fs)
    ax1.set_xlabel('Distance from Mouth of JdF (km)', fontsize=fs)

    ax1.tick_params(labelsize=fs-2) 

    ttext = infile.replace('ic_','').replace('aa_','').replace('.p','').replace('_',' ').title()
    ax1.set_title('IC: ' + ttext, fontsize=fs)
    
    ax1.text(.05,.8, 'Day = %d' % (tt),
        fontweight='bold',fontsize=fs, transform=ax1.transAxes)
    
    ax1.text(.95,.5, 'Dashed = Upper Layer',
        fontsize=fs, horizontalalignment='right', transform=ax1.transAxes)
    ax1.text(.95,.4, 'Solid = Lower Layer',
        fontsize=fs, horizontalalignment='right', transform=ax1.transAxes)
    
    if testing:
        plt.show()
    else:
        nouts = ('0000' + str(ii))[-4:]
        plt.savefig(outdir + 'plot_' + nouts + '.png')
        plt.close()
    
    ii += 1

if not testing:
    ff_str = ("ffmpeg -r 8 -i " + 
        outdir+"plot_%04d.png -vcodec libx264 -pix_fmt yuv420p -crf 25 "
        +outdir+"movie.mp4")
    os.system(ff_str)

