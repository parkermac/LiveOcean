"""
Varibles and functions used by the "flux" code.
"""
import numpy as np
from datetime import datetime, timedelta
import zfun # path provided by calling code
import pandas as pd
import Lfun
import pickle

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
# associated with lines like QQp[QQ<=0] = np.nan

def get_fluxes(indir, sect_name, in_sign=1):
    # form time series of net 2-layer transports into (+) and out of (-) the volume
    bulk = pickle.load(open(indir + 'bulk/' + sect_name + '.p', 'rb'))
    QQ = bulk['QQ']
    SS = bulk['SS']
    ot = bulk['ot']
    dt2 = []
    for tt in ot:
        dt2.append(Lfun.modtime_to_datetime(tt))
    # separate inflowing and outflowing transports
    QQp = QQ.copy()
    QQp[QQ<=0] = np.nan
    QQm = QQ.copy()
    QQm[QQ>=0] = np.nan
    # form two-layer versions of Q and S
    if in_sign == 1:
        Qin = np.nansum(QQp, axis=1)
        QSin = np.nansum(QQp*SS, axis=1)
        Qout = np.nansum(QQm, axis=1)
        QSout = np.nansum(QQm*SS, axis=1)
    elif in_sign == -1:
        Qin = -np.nansum(QQm, axis=1)
        QSin = -np.nansum(QQm*SS, axis=1)
        Qout = -np.nansum(QQp, axis=1)
        QSout = -np.nansum(QQp*SS, axis=1)
    Sin = QSin/Qin
    Sout = QSout/Qout
    fnet = bulk['fnet_lp'] # net tidal energy flux
    qabs = bulk['qabs_lp'] # low pass of absolute value of net transport
    tef_df = pd.DataFrame(index=dt2)
    tef_df['Qin']=Qin
    tef_df['Qout']=Qout
    tef_df['QSin']=QSin
    tef_df['QSout']=QSout
    tef_df['Sin']=Sin
    tef_df['Sout']=Sout
    tef_df['Ftide'] = fnet * in_sign
    tef_df['Qtide'] = qabs
    return tef_df
    
def get_budgets(sv_lp_df, v_lp_df, riv_df, tef_df_dict, seg_list):
    # volume budget
    vol_df = pd.DataFrame(0, index=sv_lp_df.index, columns=['Qin','-Qout', 'Qtide'])
    for sect_name in tef_df_dict.keys():
        df = tef_df_dict[sect_name]
        vol_df['Qin'] = vol_df['Qin'] + df['Qin']
        vol_df['-Qout'] = vol_df['-Qout'] - df['Qout']
        vol_df['Qtide'] = vol_df['Qtide'] + df['Qtide']
    vol_df['Qr'] = riv_df.sum(axis=1)
    v = v_lp_df.sum(axis=1).to_numpy()
    vol_df.loc[:,'V'] = v
    vol_df.loc[1:-1, 'dV_dt'] = (v[2:] - v[:-2]) / (2*86400)
    vol_df['Error'] = vol_df['dV_dt'] - vol_df.loc[:,'Qin'] + vol_df.loc[:,'-Qout'] - vol_df.loc[:,'Qr']
    vol_rel_err = vol_df['Error'].mean()/vol_df['Qr'].mean()

    # salt budget
    salt_df = pd.DataFrame(0, index=sv_lp_df.index, columns=['QSin','-QSout','Ftide'])
    for sect_name in tef_df_dict.keys():
        df = tef_df_dict[sect_name]
        salt_df['QSin'] = salt_df['QSin'] + df['QSin']
        salt_df['-QSout'] = salt_df['-QSout'] - df['QSout']
        salt_df['Ftide'] = salt_df['Ftide'] + df['Ftide']
    sn = sv_lp_df[seg_list].sum(axis=1).values
    salt_df.loc[1:-1, 'dSnet_dt'] = (sn[2:] - sn[:-2]) / (2*86400)
    salt_df['Smean'] = sv_lp_df.sum(axis=1)/vol_df['V']
    salt_df['Error'] = salt_df['dSnet_dt'] - salt_df['QSin'] + salt_df['-QSout']
    salt_rel_err = salt_df['Error'].mean()/salt_df['QSin'].mean()
    # add a few more columns to plot in a different way
    salt_df['Qe'] = (vol_df['Qin'] + vol_df['-Qout'])/2
    salt_df['Qnet'] = (vol_df['-Qout'] - vol_df['Qin'])
    salt_df['Sin'] = salt_df['QSin']/vol_df['Qin']
    salt_df['Sout'] = salt_df['-QSout']/vol_df['-Qout']
    salt_df['DS'] = salt_df['Sin'] - salt_df['Sout']
    salt_df['Sbar'] = (salt_df['QSin']/vol_df['Qin'] + salt_df['-QSout']/vol_df['-Qout'])/2
    salt_df['QeDS'] = salt_df['Qe'] * salt_df['DS']
    salt_df['-QrSbar'] = -salt_df['Qnet'] * salt_df['Sbar']
    salt_rel_err_qe = salt_df['Error'].mean()/salt_df['QeDS'].mean()

    # make sure everything is numeric
    for cn in vol_df.columns:
        vol_df[cn] = pd.to_numeric(vol_df[cn])
    for cn in salt_df.columns:
        salt_df[cn] = pd.to_numeric(salt_df[cn])
        
    return vol_df, salt_df, vol_rel_err, salt_rel_err, salt_rel_err_qe

# desired time ranges, the "seasons"
def get_dtr(year):
    dtr = {}
    dtr['full'] = (datetime(year,1,1,12,0,0), datetime(year,12,31,12,0,0))
    dtr['winter'] = (datetime(year,1,1,12,0,0), datetime(year,3,31,12,0,0)) # JFM
    dtr['spring'] = (datetime(year,4,1,12,0,0), datetime(year,6,30,12,0,0)) # AMJ
    dtr['summer'] = (datetime(year,7,1,12,0,0), datetime(year,9,30,12,0,0)) # JAS
    dtr['fall'] = (datetime(year,10,1,12,0,0), datetime(year,12,31,12,0,0)) # OMD
    return dtr

# Lists of 2-layer segments to use for "initial condition" experiments in the flux_engine.
# The keys should match up with "src" values in flux_engine.py.
ic_seg2_dict = {'IC_HoodCanalInner': ['H'+str(n)+'_s' for n in range(3,9)],
                                  #+ ['H'+str(n)+'_f' for n in range(3,9)],
            'IC_HoodCanal': ['H'+str(n)+'_s' for n in range(1,9)],
            'IC_SouthSound': ['S'+str(n)+'_s' for n in range(1,5)],
            'IC_SoG': ['G'+str(n)+'_s' for n in range(1,7)],
                            }
                            
season_list = list(get_dtr(2017).keys())

# create Series of two-layer volumes
# this is the one place to specify the ratio volume in the "salty" and "fresh" layers
def get_V(v_df):
    V = pd.Series()
    for seg_name in v_df.index:
        V[seg_name+'_s'] = 0.8 * v_df.loc[seg_name,'volume m3']
        V[seg_name+'_f'] = 0.2 * v_df.loc[seg_name,'volume m3']
    return V

# segment definitions, assembled by looking at the figures
# created by plot_thalweg_mean.py
segs = {
        'J1':{'S':[], 'N':[], 'W':['jdf1'], 'E':['jdf2'], 'R':['sanjuan', 'hoko']},
        'J2':{'S':[], 'N':[], 'W':['jdf2'], 'E':['jdf3'], 'R':[]},
        'J3':{'S':[], 'N':[], 'W':['jdf3'], 'E':['jdf4'], 'R':['elwha']},
        'J4':{'S':[], 'N':['sji1'], 'W':['jdf4'], 'E':['ai1','dp'], 'R':['dungeness']},
        
        'G1':{'S':['sji1'], 'N':['sji2'], 'W':[], 'E':[], 'R':['samish']},
        'G2':{'S':['sji2'], 'N':['sog1'], 'W':[], 'E':[], 'R':['nooksack', 'cowichan']},
        'G3':{'S':['sog1'], 'N':['sog2'], 'W':[], 'E':[], 'R':['nanaimo', 'fraser']},
        'G4':{'S':['sog2'], 'N':[], 'W':['sog3'], 'E':[], 'R':['clowhom', 'squamish']},
        'G5':{'S':[], 'N':['sog4'], 'W':[], 'E':['sog3'], 'R':['englishman', 'tsolum', 'oyster']},
        'G6':{'S':['sog4'], 'N':['sog5'], 'W':[], 'E':[], 'R':[]},
        
        'A1':{'S':['ai2'], 'N':[], 'W':['ai1'], 'E':[], 'R':[]},
        'A2':{'S':['ai3'], 'N':['ai2'], 'W':[], 'E':[], 'R':[]},
        'A3':{'S':['hc1'], 'N':['ai3'], 'W':[], 'E':['ai4'], 'R':[]},
        
        'M1':{'S':['mb1'], 'N':['wb1'], 'W':['ai4'], 'E':[], 'R':[]},
        'M2':{'S':['mb2'], 'N':['mb1'], 'W':[], 'E':[], 'R':[]},
        'M3':{'S':['mb3'], 'N':['mb2'], 'W':[], 'E':[], 'R':['green', 'cedar']},
        'M4':{'S':['mb4'], 'N':['mb3'], 'W':[], 'E':[], 'R':[]},
        'M5':{'S':['mb5'], 'N':['mb4'], 'W':[], 'E':[], 'R':[]},
        'M6':{'S':['tn1'], 'N':['mb5'], 'W':[], 'E':[], 'R':['puyallup']},
        
        'T1':{'S':['tn2'], 'N':['tn1'], 'W':[], 'E':[], 'R':[]},
        'T2':{'S':['tn3'], 'N':['tn2'], 'W':[], 'E':[], 'R':[]},
        
        'S1':{'S':[], 'N':['tn3'], 'W':['ss1'], 'E':[], 'R':[]},
        'S2':{'S':[], 'N':[], 'W':['ss2'], 'E':['ss1'], 'R':['nisqually']},
        'S3':{'S':[], 'N':[], 'W':['ss3'], 'E':['ss2'], 'R':[]},
        'S4':{'S':[], 'N':[], 'W':[], 'E':['ss3'], 'R':['deschutes']},
        
        'W1':{'S':['wb1'], 'N':['wb2'], 'W':[], 'E':[], 'R':['snohomish']},
        'W2':{'S':['wb2'], 'N':['wb3'], 'W':[], 'E':[], 'R':['stillaguamish']},
        'W3':{'S':['wb3'], 'N':[], 'W':[], 'E':['wb4'], 'R':[]},
        'W4':{'S':[], 'N':[], 'W':['wb4', 'dp'], 'E':[], 'R':['skagit']},
        
        'H1':{'S':['hc2'], 'N':['hc1'], 'W':[], 'E':[], 'R':[]},
        'H2':{'S':[], 'N':['hc2'], 'W':['hc3'], 'E':[], 'R':[]},
        'H3':{'S':['hc4'], 'N':[], 'W':[], 'E':['hc3'], 'R':['duckabush', 'dosewallips']},
        'H4':{'S':['hc5'], 'N':['hc4'], 'W':[], 'E':[], 'R':['hamma']},
        'H5':{'S':['hc6'], 'N':['hc5'], 'W':[], 'E':[], 'R':[]},
        'H6':{'S':[], 'N':['hc6'], 'W':[], 'E':['hc7'], 'R':['skokomish']},
        'H7':{'S':[], 'N':[], 'W':['hc7'], 'E':['hc8'], 'R':[]},
        'H8':{'S':[], 'N':[], 'W':['hc8'], 'E':[], 'R':[]},
        
        #'##':{'S':[], 'N':[], 'W':[], 'E':[], 'R':[]},
        }
        
# make lists of the various segment sequences (used below)
ssJ = ['J'+str(s) for s in range(1,5)]
ssM = ['M'+str(s) for s in range(1,7)]
ssA = ['A'+str(s) for s in range(1,4)]
ssT = ['T'+str(s) for s in range(1,3)]
ssS = ['S'+str(s) for s in range(1,5)]
ssG = ['G'+str(s) for s in range(1,7)]
ssW = ['W'+str(s) for s in range(1,5)]
ssH = ['H'+str(s) for s in range(1,9)]

# This list is the same as the keys for all the dicts below.
# we make it to have a fixed order for processing things, since
# the order of the keys of a dict may not be fixed.
channel_list = ['Juan de Fuca to Strait of Georgia',
            'Admiralty Inlet to South Sound',
            'Hood Canal',
            'Whidbey Basin']

# also cue up a line for the target salinities from the TEF sections
channel_dict = {'Juan de Fuca to Strait of Georgia':['jdf1','jdf2','jdf3','jdf4',
                'sji1', 'sji2', 'sog1','sog2','sog3','sog4','sog5'],
            'Admiralty Inlet to South Sound': ['ai1', 'ai2', 'ai3','ai4',
                'mb1','mb2','mb3','mb4','mb5',
                'tn1','tn2','tn3',
                'ss1','ss2','ss3'],
            'Hood Canal':['hc1','hc2','hc3','hc4','hc5','hc6','hc7','hc8'],
            'Whidbey Basin':['wb1','wb2','wb3','wb4','dp']}

long_channel_dict = {'Juan de Fuca to Strait of Georgia':['jdf1','jdf2','jdf3','jdf4',
                'sji1', 'sji2', 'sog1','sog2','sog3','sog4','sog5'],
            'Admiralty Inlet to South Sound': ['jdf4', 'ai1', 'ai2', 'ai3','ai4',
                'mb1','mb2','mb3','mb4','mb5',
                'tn1','tn2','tn3',
                'ss1','ss2','ss3'],
            'Hood Canal':['ai3', 'hc1','hc2','hc3','hc4','hc5','hc6','hc7','hc8'],
            'Whidbey Basin':['ai4', 'wb1','wb2','wb3','wb4','dp']}
                
seg_dict = {'Juan de Fuca to Strait of Georgia': ssJ + ssG,
            'Admiralty Inlet to South Sound': ['J4'] + ssA + ssM + ssT + ssS,
            'Hood Canal': ['A3'] + ssH,
            'Whidbey Basin': ['M1'] + ssW}
            
# same as seg_dict, but without the connections to adjoining channels
short_seg_dict = {'Juan de Fuca to Strait of Georgia': ssJ + ssG,
            'Admiralty Inlet to South Sound': ssA + ssM + ssT + ssS,
            'Hood Canal': ssH,
            'Whidbey Basin': ssW}
            
# colors to associate with each channel (the keys in channel_ and seg_dict)
clist = ['blue', 'red', 'olive', 'orange']

def make_dist(x,y):
    NS = len(x)
    xs = np.zeros(NS)
    ys = np.zeros(NS)
    xs, ys = zfun.ll2xy(x, y, x[0], y[0])
    dx = np.diff(xs)
    dy = np.diff(ys)
    dd = (dx**2 + dy**2)**.5 # not clear why np.sqrt throws an error
    dist = np.zeros(NS)
    dist[1:] = np.cumsum(dd/1000) # convert m to km
    return dist
        
def update_mm(ji, mm, this_ji_list, full_ji_list, next_ji_list):
    if mm[ji] == True:
        this_ji_list.append(ji)
        full_ji_list.append(ji)
    for ji in this_ji_list:
        mm[ji] = False
    keep_looking = True
    counter = 0
    while len(this_ji_list) > 0:
        #print('iteration ' + str(counter))
        for ji in this_ji_list:
            JI = (ji[0]+1, ji[1]) # North
            if mm[JI] == True:
                next_ji_list.append(JI)
                mm[JI] = False
            else:
                pass
            JI = (ji[0], ji[1]+1) # East
            if mm[JI] == True:
                next_ji_list.append(JI)
                mm[JI] = False
            else:
                pass
            JI = (ji[0]-1, ji[1]) # South
            if mm[JI] == True:
                next_ji_list.append(JI)
                mm[JI] = False
            else:
                pass
            JI = (ji[0], ji[1]-1) # West
            if mm[JI] == True:
                next_ji_list.append(JI)
                mm[JI] = False
            else:
                pass
        for ji in next_ji_list:
            full_ji_list.append(ji)
        this_ji_list = next_ji_list.copy()
        next_ji_list = []
        counter += 1
    return mm, this_ji_list, full_ji_list, next_ji_list