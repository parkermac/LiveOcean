"""
Varibles and functions used by the "flux" code.
"""
import numpy as np
import zfun # path provided by calling code

# segment definitions, assembled by looking at the figure
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

# also cue up a line for the target salinities from the TEF sections
channel_dict = {'Juan de Fuca to Strait of Georgia':['jdf1','jdf2','jdf3','jdf4','sji1', 'sji2', 'sog1','sog2','sog3','sog4','sog5'],
            'Admiralty Inlet to South Sound': ['ai1', 'ai2', 'ai3','ai4',
                'mb1','mb2','mb3','mb4','mb5',
                'tn1','tn2','tn3',
                'ss1','ss2','ss3'],
            'Hood Canal':['hc1','hc2','hc3','hc4','hc5','hc6','hc7','hc8'],
            'Whidbey Basin':['wb1','wb2','wb3','wb4','dp']}
                
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
clist = ['purple', 'orange', 'green', 'blue']

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