"""
The river class code.
"""

import pandas as pd
import urllib
import xml.etree.ElementTree as ET

class River:
       
    def __init__(self, Ldir):
        # initialize all expected fields
        self.att_list = ['name', 'long_name', 'usgs_code',
            'nws_code', 'scale_gage', 'scale_factor',
            'got_data', 'memo', 'flow_units']
        for att in self.att_list:
            setattr(self, att, '')
        self.Ldir = Ldir
        # default values
        self.Q = []
        self.T = []
        self.qt = pd.Series(self.Q, index=self.T)
        self.scale_factor = 1.0
        self.got_data = False
        self.memo = 'no message'
    
    def name_it(self, riv_name):
        # add the name attribute
        self.name = riv_name
      
    def get_nws_info(self, riv_name):
        # This listing has all the stations with NWS Forecasts.
        fn = self.Ldir['data'] + 'rivers/USGS_NWS_Codes.csv'
        df = pd.read_csv(fn, index_col='Name')
        self.df_nws = df
        if riv_name in df.index:
            self.long_name = df.ix[riv_name, 'Station Name']
            self.usgs_code = df.ix[riv_name, 'Station Number']
            has_nws_forecast = df.ix[riv_name, 'NWS Forecast']
            if has_nws_forecast == 'YES':
                self.nws_code = df.ix[riv_name, 'NWS ID']
        else:
            pass
   
    def get_ecy_info(self, riv_name):
        # This listing, from Mohamedali et al. (2011) has a long
        # list of rivers coming into the Salish Sea, and then associates
        # each with a USGS gage and a scaling factor.
        # We cordinate this list using riv_name, and assume it is the same
        # as the lower case version of the first work in the index (Watershed Name).
        fn = self.Ldir['data'] + 'rivers/Ecology_Scale_Factors.csv'
        df = pd.read_csv(fn, index_col='Watershed Name')
        self.df_ecy = df
        for wn in df.index:
            if wn.split()[0].lower() == riv_name:
                sg = df.ix[wn]['Scale Gage'].split()[-1]
                sg = sg[sg.find('(')+1 : sg.find(')')]
                self.scale_gage = sg
                self.scale_factor = df.ix[wn]['Scale Factor']
                break
                 
    def print_info(self):
        print(50*'-')
        for att in self.att_list:
            print(att + ' = ' + str(getattr(self,att)))
        print(50*'-')
        
    def get_usgs_data(self, days):
        # This gets USGS data for a past time period specfied by
        # the tuple of datetimes "days".  If "days" is empty
        # then we get the most recent 6 days.
        
        # Set the time period to get.
        if len(days) == 2:
            time_str = ('&startDT=' + days[0].strftime('%Y-%m-%d')
                +'&endDT=' + days[1].strftime('%Y-%m-%d'))            
        else:
            # This gets the most recent 6 days (daily intervals?)
            time_str = '&period=P6D'
        gage = self.usgs_code
        if len(self.scale_gage) > 0 and self.scale_gage != 'Hybrid':
            gage = self.scale_gage
        # Form the url.
        url_str = ('http://waterservices.usgs.gov/nwis/dv/'
            + '?format=waterml,1.1&sites=' + str(gage)
            + time_str + '&parameterCd=00060')
        # Get the XML.     
        file = urllib.request.urlopen(url_str, timeout=10)
        try:
            tree = ET.parse(file)
            root = tree.getroot()
        except:
            self.memo = 'problem downloading XML'
        # Extract the data from the XML.
        try:
            flag = True
            # aa = '{http://www.cuahsi.org/waterML/1.1/}'
            rt = root.tag  
            aa = rt[rt.find('{'): rt.find('}') + 1]
            for e0 in root.findall(".//"):
                if e0.tag == aa+'value':
                    self.Q.append(float(e0.text))
                    self.T.append(pd.to_datetime(e0.get('dateTime')))
                if e0.tag == aa+'unitCode' and flag:
                    self.flow_units = e0.text
                    flag = False
            self.qt = pd.Series(self.Q, index=self.T)
            self.fix_units()
            self.qt = float(self.scale_factor) * self.qt
            # Note: this resampling of daily data just moves the timestamp to noon
            # of the day it started at.  Data is unchanged.
            self.qt = self.qt.resample('D', how='mean', label='right', loffset='-12h')
            self.got_data = True
            self.memo = 'success'
        except:
            self.memo = 'problem parsing data from XML'
        self.memo = (self.memo + ' USGS')
            
    def get_nws_data(self):
        # This gets NWS forecast data.
            
        url_str = ('http://www.nwrfc.noaa.gov/xml/xml.cgi?id=' +
                   self.nws_code +
                   '&pe=HG&dtype=b&numdays=10')
        file = urllib.request.urlopen(url_str, timeout=10)
        try:
            tree = ET.parse(file)
            root = tree.getroot()
        except:
            self.memo = 'problem downloading XML'
        
        try:
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
                                    self.Q.append(float(ee.text))
                                    if flag:
                                        self.flow_units = ee.get('units')
                                        flag = False
                                if ee.tag == aa+'dataDateTime':
                                    self.T.append(pd.to_datetime(ee.text))
            self.qt = pd.Series(self.Q, index=self.T)
            self.fix_units()
            self.qt = float(self.scale_factor) * self.qt
            self.qt = self.qt.resample('D', how='mean', label='right', loffset='-12h')
            self.got_data = True
            self.memo = 'success'
        except:
            self.memo = 'problem parsing data from XML'
        self.memo = (self.memo + ' NWS')
        
    def get_ec_data(self, days):
        # gets Environment Canada data, like for the Fraser
        try:
            indir = self.Ldir['data'] + 'rivers/data_processed/'
            s0 = pd.read_pickle(indir + self.name + '_flow_historical.p')
            s1 = pd.read_pickle(indir + self.name + '_flow_historical_2.p')
            s2 = pd.read_pickle(indir + self.name + '_flow.p')
            # trim s1 so there is no overlap
            s1 = s1[s1.index > s0.index[-1]]
            s1 = s1[s1.index < s2.index[0]]
            ss = pd.concat([s0, s1, s2])
            ss = ss[ss.index >= days[0]]
            ss = ss[ss.index <= days[1]]
            self.qt = ss
            self.got_data = True
            self.memo = 'success'
        except:
            self.memo = 'problem accessing data'
        self.memo = (self.memo + ' EC')
      
    def fix_units(self):
        # fix units
        if self.flow_units == 'kcfs':
            self.qt = self.qt*28.3168466
            self.flow_units = '$m^{3}s^{-1}$'
        elif (self.flow_units == 'cubic feet per second'
            or self.flow_units == 'ft3/s'):
            self.qt = self.qt*0.0283168466
            self.flow_units = '$m^{3}s^{-1}$'
        else:
            pass


