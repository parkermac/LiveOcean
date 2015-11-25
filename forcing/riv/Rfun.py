"""
These are functions used by the river code.
"""
    
def fix_units(qt, flow_units):
    # fix units
    if flow_units == 'kcfs':
        qt = qt*28.3168466
    elif flow_units == 'cubic feet per second' or flow_units == 'ft3/s':
        qt = qt*0.0283168466
    return qt
    
def get_nws_data(id):
    """
    This gets NWS data.    
    """
    import urllib.request as U
    import xml.etree.ElementTree as ET
    import pandas as pd
      
    # default values
    Q = []
    T = []
    qt = pd.Series(Q, index=T)
    got_nws_data = False
    memo = 'no message'
         
    url_str = ('http://www.nwrfc.noaa.gov/xml/xml.cgi?id=' + id
    + '&pe=HG&dtype=b&numdays=10')
    try:
        file = U.urlopen(url_str, timeout = 10)
        tree = ET.parse(file)
        root = tree.getroot()
    except:
        memo = 'problem downloading XML'
    
    try:
        flow_units = ''
        flag = True
        # NOTE: you find this tag by looking at any instance of e0.tag
        #aa = '{http://www.nwrfc.noaa.gov/xml/schemas/2004/03/hydromet_data}'
        rt = root.tag  
        aa = rt[rt.find('{'): rt.find('}') + 1]
        for e0 in root.findall(".//"):
            if e0.tag == aa+'observedData' or e0.tag == aa+'forecastData':
                for e in e0:
                    if e.tag == aa+'observedValue' or e.tag == aa+'forecastValue':
                        for ee in e:
                            if ee.tag == aa+'discharge':
                                Q.append(float(ee.text))
                                if flag:
                                    flow_units = ee.get('units')
                                    flag = False
                            if ee.tag == aa+'dataDateTime':
                                T.append(pd.to_datetime(ee.text))
        qt = pd.Series(Q, index=T)
        qt = fix_units(qt, flow_units)
        got_nws_data = True
        memo = 'success'
    except:
        memo = 'problem parsing data from XML'
        
    return qt, got_nws_data, memo
    
def get_usgs_data(id):
    """
    This gets USGS data.
    """
    import urllib.request as U
    import xml.etree.ElementTree as ET
    import pandas as pd
      
    # default values
    Q = []
    T = []
    qt = pd.Series(Q, index=T)
    got_usgs_data = False
    memo = 'no message'
    
    url_str = ('http://waterservices.usgs.gov/nwis/iv/' +
    '?format=waterml,1.1&sites=' + id + '&period=P6D&parameterCd=00060')
    try:
        file = U.urlopen(url_str, timeout = 10)
        tree = ET.parse(file)
        root = tree.getroot()
    except:
        memo = 'problem downloading XML'
                                          
    try:
        flow_units = ''
        flag = True
        #aa = '{http://www.cuahsi.org/waterML/1.1/}'
        rt = root.tag  
        aa = rt[rt.find('{'): rt.find('}') + 1]
        for e0 in root.findall(".//"):
            if e0.tag == aa+'value':
                Q.append(float(e0.text))
                T.append(pd.to_datetime(e0.get('dateTime')))
            if e0.tag == aa+'unitCode' and flag:
                flow_units = e0.text
                flag = False
        qt = pd.Series(Q, index=T)
        qt = fix_units(qt, flow_units)
        got_usgs_data = True
        memo = 'success'
    except:
        memo = 'problem parsing data from XML'
    
    return qt, got_usgs_data, memo
    
def get_usgs_data_past(id, datetime_start, datetime_end):
    """
    This gets USGS data for some past time period.
    """
    import urllib.request as U
    import xml.etree.ElementTree as ET
    import pandas as pd
      
    # default values
    Q = []
    T = []
    qt = pd.Series(Q, index=T)
    got_usgs_data = False
    memo = 'no message'
    
    # this gives daily values from a start date to an end date
    # note that parameterCd=00010 is temperature degC (separate parameters with a comma)
    url_str = ('http://waterservices.usgs.gov/nwis/dv/'
        + '?format=waterml,1.1&sites=' + id + '&startDT='
        + datetime_start.strftime('%Y-%m-%d')
        +'&endDT=' + datetime_end.strftime('%Y-%m-%d') + '&parameterCd=00060')
    
    try:
        file = U.urlopen(url_str, timeout = 10)
        tree = ET.parse(file)
        root = tree.getroot()
    except:
        memo = 'problem downloading XML'
                                          
    try:
        flow_units = ''
        flag = True
        #aa = '{http://www.cuahsi.org/waterML/1.1/}'
        rt = root.tag  
        aa = rt[rt.find('{'): rt.find('}') + 1]
        for e0 in root.findall(".//"):
            if e0.tag == aa+'value':
                Q.append(float(e0.text))
                T.append(pd.to_datetime(e0.get('dateTime')))
            if e0.tag == aa+'unitCode' and flag:
                flow_units = e0.text
                flag = False
        qt = pd.Series(Q, index=T)
        qt = fix_units(qt, flow_units)
        got_usgs_data = True
        memo = 'success'
    except:
        memo = 'problem parsing data from XML'
    
    return qt, got_usgs_data, memo


def get_river_code_nws(river_name, Ldir):
    """
    This goes from a river name to an NWS number.
    """
    
    import pandas as pd

    # get the table of names, gages, and scaling factors from Mohamedali
    inframe = pd.read_csv(Ldir['LO'] + 'forcing/riv/' +
        'Files_USGS/USGS_good_plus.csv', index_col='Name')

    try:
        has_nws_forecast = inframe.ix[river_name, 'NWS Forecast'] # YES or NO
        this_code_nws = inframe.ix[river_name, 'NWS ID']
        this_code_usgs = str(inframe.ix[river_name, 'Station Number']) 
    except:
        has_nws_forecast = 'NO'
        this_code_nws = 'MISSING'
        this_code_usgs = 'MISSING'
   
    return this_code_nws, has_nws_forecast, this_code_usgs
    
def get_river_code_ecology(river_name, Ldir):
    """
    This goes from a river name
    to a USGS number and a scale factor, based on the
    spreadsheet by Mohamedali et al. at Ecology.
    """
        
    import pandas as pd

    # get the table of names, gages, and scaling factors from Mohamedali
    inframe = pd.read_csv(Ldir['LO'] + 'forcing/riv/' +
        'Files_Ecology/Ecology_Scale_Factors.csv')
    name = inframe['Watershed Name'].values
    gage = inframe['Scale Gage'].values
    # fac = inframe['Scale Factor'].values

    # pull the station code out
    # issue: most are USGS numbers, but some are 'Hybrid' or like '08HA070' (Canadian)
    # but all are still string items in a list
    gg = []
    for g in gage: # find the item in parentheses
        i1 = g.find('(')
        i2 = g.find(')')
        gg.append(g[i1+1:i2])

    # simplify the names   
    rn = []
    for n in name:
        n = n.strip()
        i1 = n.find(' ')
        if i1 != -1:
            n = n[:i1]
        i2 = n.find('_')
        if i2 != -1:
            n = n[:i2]
        rn.append(n.lower())
        
    # add these to the frame
    inframe['Code'] = gg

    # and index by River Name
    if2 = inframe.set_index([rn])
       
    try:
        this_code = if2.ix[river_name, 'Code']
        this_scale = if2.ix[river_name, 'Scale Factor']
    except:
        this_code = 'MISSING'
        this_scale = 1
    
    return this_code, this_scale
    
def get_rivers_dataframe(Ldir):
    """
    Gets a dataframe with all the needed info for all the rivers.
    """
    import pandas as pd
       
    # get the list of rivers that we need for this run
    rnames_frame = pd.read_csv(Ldir['run'] + 'rname_list.txt', header=None,
        names=['River Name'])
    rnames_full = rnames_frame['River Name'].values
    rnames_full = rnames_full.tolist()
    
    # swap some names
    try:
        rnames_full[rnames_full.index('duwamish')] = 'green'
        rnames_full[rnames_full.index('hammahamma')] = 'hamma'
    except:
        print('Problem swapping river names')
    
    # DEBUGGING
    if True:
        # use a list of all the rivers
        rnames = rnames_full
    else:
        # or just some
        rnames = ['skagit',  'columbia', 'fraser']
    
    # initialize the lists
    usgs_id = []
    nws_id = []
    can_id = []
    scale_factor = []
    
    print_stuff = False
    if print_stuff: # helpful for debugging, first print a header line
        print('%15s %10s %10s %10s %10s %10s %10s' % ('name', 'USGS_ecy', 'ECY_sf',
        'NWS_code', 'Has_NWS', 'USGS_nws', 'code_usgs'))
    
    # step through the river names    
    for rn in rnames_full:
        
        # this function reads a file that has scale factors linked to USGS numbers (right?)
        code_ecology, scale_factor_ecology = get_river_code_ecology(rn, Ldir)
        
        # this function reads a file that gives NWS (forecast) codes for some USGS gauged rivers
        code_nws, has_nws_forecast, code_nws_usgs = get_river_code_nws(rn, Ldir)
        
        # set code_usgs  
        if code_ecology == code_nws_usgs:
            code_usgs = code_nws_usgs
        elif (code_ecology == 'MISSING' or code_ecology == 'Hybrid') and code_nws_usgs != 'MISSING':
            code_usgs = code_nws_usgs
        elif code_nws_usgs == 'MISSING' and (code_ecology != 'MISSING' and code_ecology != 'Hybrid'):
            code_usgs = code_ecology
        else:
            code_usgs = 'GACK'
        
        # intervene by hand for skagit_south   
        if rn == 'skagit_south':
            code_usgs = '12200500'
            scale_factor_ecology = 0
        
        # set code_can (Canadian)    
        if rn == 'fraser':
            code_can = code_ecology
            code_usgs = 'MISSING'
        else:
            code_can = 'MISSING'
        
        # set code_nws 
        if has_nws_forecast == 'YES':
            code_nws = code_nws
        elif has_nws_forecast == 'NO':
            code_nws = 'MISSING'
        else:
            code_nws = 'GACK'
        
        # then add items to the lists    
        usgs_id.append(code_usgs)
        nws_id.append(code_nws)
        can_id.append(code_can)
        scale_factor.append(scale_factor_ecology)
        
        if print_stuff: # helpful for debugging
            print('%15s %10s %10.2f %10s %10s %10s %10s' % (rn, code_ecology,
                scale_factor_ecology, code_nws, has_nws_forecast,
                code_nws_usgs, code_usgs))
    
    # put the results in a DataFrame
    rdf = pd.DataFrame(index=rnames_full)
    rdf.index.name = 'River_Name'
    rdf['NWS_ID'] = nws_id
    rdf['USGS_ID'] = usgs_id
    rdf['CAN_ID'] = can_id
    rdf['Scale_Factor'] = scale_factor
    # tidy!
    
    return rdf, rnames
    
def get_rivers_data(rdf, rnames, Info, Ldir):
    """
    Go out and get the data for all.
    """
    import numpy as np
    from pandas import DataFrame
    import pandas as pd
    from datetime import datetime
    
    # next we go through the rivers in rdf and get the data for each
    riv_dict = {}
    good_riv_list = []
    bad_riv_list = []
    
    for rn in rnames: #rdf.index:
        # set up some things
        got_nws_data = False
        got_usgs_data = False
        
        # get the ID's and scale factor for this river
        aa = rdf.ix[rn]
        
        print('** working on ' + rn)
        
        if Info['run_type'] == 'forecast':
            # get the NWS Forecast
            if aa['NWS_ID'] != 'MISSING':
                id = aa['NWS_ID']
                print(' >> trying NWS: ' + id)
                qt, got_nws_data, memo = get_nws_data(id)
                print('  -----' + memo + '-------')
                # and if that fails try the USGS  recent observations
                if got_nws_data == False:
                    id = aa['USGS_ID']
                    print('  >> trying USGS: ' + id)
                    qt, got_usgs_data, memo = get_usgs_data(id)
                    print('  -----' + memo + '-------')
        elif Info['run_type'] == 'backfill':
            datetime_start = datetime.strptime(Info['datestring_start'],'%Y.%m.%d')
            datetime_end = datetime.strptime(Info['datestring_end'],'%Y.%m.%d')
            id = aa['USGS_ID']
            # really we should have a trap door here if id == 'MISSING'
            # but it is handled semi-gracefully anyway
            print('  >> trying USGS backfill: ' + id)
            qt, got_usgs_data, memo = get_usgs_data_past(id,
                datetime_start, datetime_end)
            print('  -----' + memo + '-------')
        # then do some processing and pack them in a dict
        if got_nws_data | got_usgs_data:  
            # fix time zone (e.g. USGS reports in local time)    
            # qt = qt.tz_convert('UTC')
            # 3/21/2014 Why doens't this work anymore - as of move to Canopy
            # block average to daily values
            qt = qt.resample('D', how='mean', label='right', loffset='-12h')
            # fix the scaling
            qt = qt*aa['Scale_Factor']
            
            # pack into an item in a dict
            riv_dict[rn] = qt
            good_riv_list.append(rn)
        else:
            bad_riv_list.append(rn)
            
    # convert riv_dict to a DataFrame (plot using riv_df.plot() - very cool!)
    riv_df = DataFrame(riv_dict)
    
    # make sure all values are filled (assumes persistence is OK)
    riv_df = riv_df.fillna(method='bfill')
    riv_df = riv_df.fillna(method='ffill')
    
    # add columns for the missing rivers (will fill with climatology)
    # there has to be a cleaner way to do this!
    missing_riv_list = [val for val in rnames if val not in good_riv_list]
    for rn in missing_riv_list:
        riv_df[rn] = np.nan
    # riv_df should now have columns for ALL the rivers, many with no data
    
    # load the climatology
    # for the flow
    clim_fn = (Ldir['LO'] + 'forcing/riv/' +
        'river_climatology/Qclim.csv')
    clim_df = pd.read_csv(clim_fn, index_col='yearday')
    # and for temperature
    tclim_fn = (Ldir['LO'] + 'forcing/riv/' +
        'river_climatology/Tclim.csv')
    tclim_df = pd.read_csv(tclim_fn, index_col='yearday')
    # note clim_df.plot() makes a nice plot of this!
    
    # use climatology for missing values
    riv_df_new = riv_df.copy()
    for rn in riv_df:
        this_riv = riv_df[rn]
        isgood = this_riv.isnull()
        for t in this_riv.index:
            if isgood[t]:
                yd = t.dayofyear
                cval = clim_df.ix[yd, rn] * rdf.ix[rn, 'Scale_Factor']
                riv_df_new.ix[t, rn] = cval
                
    # make temperature purely from climatology
    triv_df = DataFrame(np.nan, index=riv_df_new.index,
        columns=riv_df_new.columns)
    triv_df_new = triv_df.copy()
    for rn in triv_df:
        this_triv = triv_df[rn]
        isgood = this_triv.isnull()
        for t in this_triv.index:
            if isgood[t]:
                yd = t.dayofyear
                cval = tclim_df.ix[yd, rn]
                triv_df_new.ix[t, rn] = cval
    
    # make sure that the columns are in the correct order
    riv_df_new = riv_df_new.reindex(columns=rnames)
    triv_df_new = triv_df_new.reindex(columns=rnames)
    
    # it would be nice to check here that the result is good (no NaN's)
    
    # create a DataFrame that uses seconds since 1/1/1970 as its index
    # add a column of seconds since 1970/01/01
    tt = riv_df_new.index   
    TT = []
    for ttt in tt:    
        TT.append(ttt.value)
    T = np.array(TT)
    Tsec = T/1e9
    
    riv_df_for_matlab = DataFrame(riv_df_new.values, index=Tsec,
        columns=riv_df_new.columns)
    riv_df_for_matlab.index.name = 'Tsec'
    
    triv_df_for_matlab = DataFrame(triv_df_new.values, index=Tsec,
        columns=triv_df_new.columns)
    triv_df_for_matlab.index.name = 'Tsec'
    
    return riv_df_for_matlab, triv_df_for_matlab






