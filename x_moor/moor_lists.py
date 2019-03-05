"""
Module to create dicts for multiple (or single) mooring extractions.
"""

def get_sta_dict(job_name):
    
    # Initialize lists of varibles to get.  If you leave them empty
    # then the calling code will get everything.
    v2_list = [] # 2-D variables (like 'zeta')
    v3_list_rho = [] # 3-D variables on the vertical rho grid (like 'salt')
    v3_list_w = [] # 3-D varibles on the vertical w grid (like 'w')
    
    # specific job definitions
    
    if job_name == 'willapa_bc': # Willapa Bay Center PCSGA Mooring
        sta_dict = {
            'wbc': (-123.9516, 46.6290)
            }
            
    elif job_name == 'kd_array': # Kathy Donohue, Slow Slip related 2018_09
        sta_dict = {
            'kda00': (-125.155750, 44.447250),
            'kda01': (-125.206174, 44.447359),
            'kda02': (-125.256597, 44.447469),
            'kda03': (-125.307021, 44.447578),
            'kda04': (-125.357446, 44.447688),
            'kda05': (-125.407870, 44.447797),
            'kda06': (-125.458295, 44.447906),
            'kda07': (-125.155750, 44.402250),
            'kda08': (-125.206134, 44.401207),
            'kda09': (-125.256516, 44.400165),
            'kda10': (-125.306896, 44.399123),
            'kda11': (-125.357275, 44.398081),
            'kda12': (-125.407651, 44.397039),
            'kda13': (-125.458026, 44.395997),
            'kda14': (-125.155750, 44.357249),
            'kda15': (-125.206095, 44.355605),
            'kda16': (-125.256436, 44.353961),
            'kda17': (-125.306775, 44.352317),
            'kda18': (-125.357111, 44.350672),
            'kda19': (-125.407445, 44.349028),
            'kda20': (-125.457775, 44.347384),
            'kda21': (-125.085750, 44.177247),
            'kda22': (-125.135941, 44.175703),
            'kda23': (-125.186129, 44.174159),
            'kda24': (-125.236314, 44.172614),
            'kda25': (-125.286497, 44.171070),
            'kda26': (-125.336677, 44.169525),
            'kda27': (-125.386855, 44.167981)
            }
        v2_list = ['zeta', 'Pair']
        v3_list_rho = ['rho', 'u', 'v']
        v3_list_w = []
        
    elif job_name == 'pelletier_1':
        sta_dict = {
            'NANOOS_ORCA_Carr_Inlet': (-122.730000,47.2800),
            'NANOOS_ORCA_Hoodsport': (-123.112600,47.4218),
            'NANOOS_ORCA_Point_Wells': (-122.397200,47.7612),
            'NANOOS_ORCA_Hansville': (-122.627000,47.9073),
            'NANOOS_ORCA_Twanoh': (-123.008300,47.3750),
            'NANOOS_ORCA_Dabob_Bay': (-122.802900,47.8034),
            'UW_Bellingham_Bay_Buoy': (-122.576500,48.7237),
            'NANOOS_ChaBa_Buoy': (-124.950000,47.9700),
            'NANOOS_ORCA_Duckabush': (-123.000000,47.5533),
            'P12': (-123.108000,47.4253),
            'P22': (-123.019000,48.2717),
            'P28': (-122.454000,47.7034),
            'P38': (-122.708000,47.2766),
            'P4': (-122.553000,48.2422),
            'P402': (-123.023000,47.3567),
            'P8': (-122.605000,47.8967),
        }
        
    elif job_name == 'pelletier_2':
        sta_dict = {
            'P381': (-124.949, 47.969),
            'P136': (-123.483, 48.22333)
        }
        
    elif job_name == 'erika_1':
        sta_dict = {
        'heinBankBuoy': (-123.165000, 48.334000),
        'portWilliamsBuoy': (-122.406100, 47.537200),
        'offshore': (-126.426364, 47.595087),
        'jdfStraitEntrance': (-124.788896, 48.473260),
        'hoodCanalSouth': (-122.762032, 47.716531),
        'hoodCanalNorth': (-122.603685, 47.871478),
        'hoodsportBuoy': (-123.112600, 47.421800),
        'deceptionPassEast': (-122.599176, 48.414963),
        'deceptionPassWest': (-122.688143, 48.400888),
        'bellinghamBayBuoySite': (-122.576500, 48.723700),
        'saratogaPassage': (-122.552402, 48.156914),
        'admiraltyNorth': (-122.851869, 48.196848),
        'admiraltySouth': (-122.656958, 48.094407),
        'bellinghamBaySouth': (-122.541997, 48.639757),
        'seattle': (-122.455031, 47.559993),
        'jdfStraitMiddle': (-123.803274, 48.246187)
        }
        
    elif job_name == 'comt3_2014_offshore':
        sta_dict = {
            'CA015': (-124.756833,48.166300),
            'CA042': (-124.823367,48.166017),
            'CE015': (-124.349050,47.354167),
            'CE042': (-124.488733,47.353133),
            'KL015': (-124.428400,47.600833),
            'KL027': (-124.497067,47.594567),
            'MB015': (-124.681883,48.325433),
            'MB042': (-124.735383,48.323967),
            'TH015': (-124.618650,47.875500),
            'TH042': (-124.733417,47.876150),
            'CHABA': (-124.949233,47.965950),
            'NEMO': (-124.948300,47.963333),
            'NH10_CARBON': (-124.300000,44.650000),
            'NH20_CARBON': (-124.500000,44.650000),
            'CE06ISSM': (-124.272000,47.133600),
            'CE07SHSM': (-124.566000,46.985900),
            'CE09OSSM': (-124.972000,46.850800),
            'CE01ISSM': (-124.095000,44.659800),
            'CE02SHSM': (-124.304000,44.639300),
            'CE04OSSM': (-124.956000,44.381100),
            'CE02SHBP': (-124.306000,44.637100),
            'CE04OSBP': (-124.954000,44.369500),
        }
        
    return sta_dict, v2_list, v3_list_rho, v3_list_w