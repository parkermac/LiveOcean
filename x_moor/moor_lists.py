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
        
    return sta_dict, v2_list, v3_list_rho, v3_list_w