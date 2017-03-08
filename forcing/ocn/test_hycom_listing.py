import Ofun
from importlib import reload
reload(Ofun)

# # create a list of url's of the preprocessed HYCOM files for this forecast
# fn_list = Ofun.get_hycom_file_list()

import xml.etree.ElementTree as ET
import urllib.request as U

# initiate the file list
fn_list = []

# get the xml of the catalog
#xml_name = 'http://tds.hycom.org/thredds/catalog/GLBu0.08/expt_91.2/forecasts/catalog.xml'
xml_name = 'http://beta.hycom.org/thredds/catalog/GLBu0.08/expt_91.2/forecasts/catalog.xml'
try:
    xfile = U.urlopen(xml_name, timeout=30)
except:
    print('problem getting xfile')
tree = ET.parse(xfile)
xfile.close()
root = tree.getroot()

# NOTE: you find "xmlns" by looking at any instance of e0.tag (the part contained in {})
# (just type e0.tag at the command line) or by looking at the "xmlns" attribute
# in the first line of the XML listing when viewed in Chrome.
# Or, here we automate the procedure (best)...
rt = root.tag
xmlns = rt[rt.find('{'): rt.find('}') + 1]

# NOTE: in order to successfully find things in the XML you have to look at
# it first, which you do by going to the url given by "xml_name" above
# in Firefox or Chrome (Chrome is most complete)

# .// selects all subelements, on all levels beneath the current element.
# For example, .//egg selects all egg elements in the entire tree.

# get the url prefix
for e0 in root.findall('.//' + xmlns + 'service'):
    if e0.get('name') == 'ncdods':
        url_prefix = e0.get('base')

# get the remainder of the file paths and put them in a list
for e0 in root.findall('.//' + xmlns + 'dataset'):
    if e0.get('urlPath') != None:
        fn_list.append(url_prefix + e0.get('urlPath'))
        


# # get a selection of the raw list (e.g. one per day)
varf_dict, dt_list2 = Ofun.get_varf_dict(fn_list)
var_list = list(varf_dict.keys())
NT = len(dt_list2)



    
    
        