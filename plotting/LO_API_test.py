import json
import urlparse
import requests
 
# Establish the api base URL and the first call,
# 'get-info' that tells what data are available
baseUrl = "http://liveocean.azurewebsites.net/api/"
requestcall = "get-info"
 
apiReply = requests.get(urlparse.urljoin("http://liveocean.azurewebsites.net/api/", "get-info"))
print 'Call =', apiReply.url, '; reply status is', apiReply.ok, '; status code is', apiReply.status_code, '\n'
 
# convert this reply to json and print some diagnostics
rjson = apiReply.json()
print 'jsonized get-info reply is of', type(rjson), '; length', len(rjson), '; PKeys (unicode) are',
print [element['PartitionKey'] for element in rjson]
 
# Here the unicode u'xxxx' is encoded as ascii to make the PartitionKey list easier to read
pk=[]
for i in range(0,len(rjson)):
        pk.append(rjson[i]['PartitionKey'].encode('ascii','ignore'))
print ""
print 'As ascii:',pk,'\n'
 
# We move on now to the get-value api call: JSON format, contents are data (in this case salinity at 10 meters depth)
requestcall = "get-value"
data_args = {'date':'2015-08-26T01:00:00Z', 'depthMeters':'-10.0','param':'salt'}
apiReply = requests.get(urlparse.urljoin(baseUrl, requestcall), params=data_args)
 
print 'Call =', apiReply.url
print 'Reply status is', apiReply.ok
print 'Status code is', apiReply.status_code, '\n'
 
# convert response to json
jsonReply = apiReply.json()
print 'jsonized get-value reply is of', type(jsonReply), '; length', len(jsonReply)
 
# pull out the dictionary entries, 'min', 'max' and 'data'
minValue = jsonReply['min']
maxValue = jsonReply['max']
print 'minValue is of', type(minValue), 'with value', minValue
print 'maxValue is of', type(maxValue), 'with value', maxValue
print ""
print ""
data=jsonReply['data']
print 'The data entry of this dictionary is of', type(data)
 
# Determine the dimensions of the data List-of-Lists
numLists = len(data)
numElementsPerList = len(data[0])
print 'API getvalue provides', numLists, 'lists each with', numElementsPerList, 'elements','\n'
 
# Print a couple of exemplars, one from the water and one from the masked-out land
print 'element 341 of list 131, data[130][340], is',data[130][340],'which is ',type(data[130][340]),'\n'
print 'element 381 of list 174, data[173][380], is',data[173][380],'which is ',type(data[173][380]),'\n'
 