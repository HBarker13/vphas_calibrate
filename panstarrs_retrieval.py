#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python
#Retrieve panstarrs data and save
#Copied mostly from mast examples/tutorial page

import sys
import os
import time
import re
import json
from astropy.table import Table
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u


#MAST documentation here: https://mast.stsci.edu/api/v0/
# Examples https://mast.stsci.edu/api/v0/pyex.html
#Web search: http://archive.stsci.edu/panstarrs/search.php
#Tutorial: https://mast.stsci.edu/api/v0/MastApiTutorial.html


try: # Python 3.x
	from urllib.parse import quote as urlencode
	from urllib.request import urlretrieve
except ImportError:  # Python 2.x
	from urllib import pathname2url as urlencode
	from urllib import urlretrieve

try: # Python 3.x
	import http.client as httplib 
except ImportError:  # Python 2.x
	import httplib               



## [Mast Query]
def mastQuery(request):

    server='mast.stsci.edu'

    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)
    
    # opening the https connection
    conn = httplib.HTTPSConnection(server)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    # Close the https connection
    conn.close()

    return head,content





## [Json to astropy]
def mastJson2Table(jsonObj):

    dataTable = Table()

    for col,atype in [(x['name'],x['type']) for x in jsonObj['fields']]:
        if atype=="string":
            atype="str"
        if atype=="boolean":
            atype="bool"
        dataTable[col] = np.array([x.get(col,None) for x in jsonObj['data']],dtype=atype)
        
    return dataTable




## [Cone Search]   
#I emailed the help desk and got the response:
#There is no specific service to access panstarrs data.  Instead you would use the general services, Mast.Caom.Cone, Mast.Caom.Filtered, and Mast.Caom.Filtered.Position.  To perform a cone search specifically on panstarrs data, you will want the Mast.Caom.Filtered.Position services (see https://mast.stsci.edu/api/v0/_services.html#MastCaomFilteredPosition).
#
#An example query would be:
#request = {"service":"Mast.Caom.Filtered.Position",
#           "format":"json",
#           "params":{
#               "columns":"*",
#               "filters":[{"paramName":"obs_collection",  "values":["PS1"] }],
#               "position":"210.8023, 54.349, 0.0"
#           }}
#-Clara brasseur
def ConeSearch():

	request = {"service":"Mast.Caom.Filtered.Position",
		"format":"json",
	        "params":{
		"columns":"*",
	        "filters":[{"paramName":"obs_collection",  "values":["PS1"] }],
                "position": str(ra)+","+str(dec)+","+str(radius) 
           }}

	headers,outString = mastQuery(request)

	outData = json.loads(outString)

	return outData





#input ra and dec in sexadecimal
ra = raw_input('Enter pointing RA (hh:mm:ss): ')
dec = raw_input('Enter pointing Dec (+dd:mm:ss): ')

#convert to degrees
coords = SkyCoord(ra+' '+dec, frame='icrs', unit=(u.hourangle, u.deg))
ra = str(coords.ra.degree)
dec = str(coords.dec.degree)

radius = 1.2

print
print 'Searching: '
print 'RA ', ra
print 'Dec ', dec
print 'Radius ', radius
print


#create a panstarrs directory if it doesn't exists
pan_dir = os.getcwd() + '/panstarrs_data'
if not os.path.exists(pan_dir):
	print 'Creating', pan_dir
	os.makedirs(pan_dir)



pan_json = ConeSearch()

#convert the cone search into an astropy table
pan_astropy = mastJson2Table(pan_json)


#see what data products are available
for ind in range(len(pan_astropy)):

	obs = pan_astropy[ind]
	print 'Observation ', [obs[x] for x in ['dataproduct_type', 'obs_collection', 'instrument_name']]

	obsid = obs['obsid']
	productRequest = {'service':'Mast.Caom.Products',
				'params':{'obsid':obsid},
				'format':'json',
				'pagesize':10000,
				'page':1}
				
	headers,obsProductsString = mastQuery(productRequest)
	obsProducts = json.loads(obsProductsString)

	
	#get the catalgoue data. Start by making an astropy table containing just the catalogue product info
	catProdArr = [x for x in obsProducts['data'] if x.get("productType",None) == 'CATALOG']
	if len(catProdArr)==0: continue #no catalogues


	print 'Retrieving ', ind+1,'/',len(pan_astropy)
	
	catProducts = Table()
	for col,atype in [(x['name'],x['type']) for x in obsProducts['fields']]:
		if atype=="string":
        		atype="str"
    		if atype=="boolean":
        		atype="bool"
    		if atype == "int":
        		atype = "float" # array may contain nan values, and they do not exist in numpy integer arrays
    		catProducts[col] = np.array([x.get(col,None) for x in catProdArr],dtype=atype)
	
	
	#download products, this may take some time
	urls = catProducts['dataURI']
	descriptions = catProducts['description'] 
	productTypes = catProducts['dataproduct_type']
	outPaths = ["mastFiles/"+x['obs_collection']+'/'+x['obs_id']+'/'+x['productFilename'] for x in catProducts]
	zipFilename = "mastDownload"
	extension = "tar.gz"
	
	bundleRequest = {"service":"Mast.Bundle.Request",
                 "params":{"urlList":",".join(urls),
                           "filename":zipFilename,
                           "pathList":",".join(outPaths),
                           "descriptionList":list(descriptions),
                           "productTypeList":list(productTypes),
                           "extension":extension},
                 "format":"json",
                 "page":1,
                 "pagesize":100000}  

	headers,bundleString = mastQuery(bundleRequest)
	bundleInfo = json.loads(bundleString)

	savename = pan_dir+"/"+zipFilename+str(ind)+"."+extension
	print 'Saving', savename
	urlretrieve(bundleInfo['url'], savename)
	print 'Dowloaded'
	raw_input('stopped')



	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
