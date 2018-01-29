#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python
import urllib, urlparse, string, time
import os
import glob
from astropy.io import fits, ascii
from astropy import wcs
from astropy.table import vstack, unique, Table
import numpy as np
import csv
 
#Use a HTML GET to search for panstarrs data 
#Doesn't get round the search radius limits and object number limit so will have to do one ccd at a time



#skip the whole thing if the final file already exists
fin_fpath = os.getcwd()+'/panstarrs.fits'
if os.path.exists(fin_fpath):
	print 'Final file already exists'
	print fin_fpath
	import sys
	sys.exit()




#get the central coordinats of each CCD and use those coordinates in 
#a pan-starrs query. The maximum radius of a query is 30arcmin, but the
#object limit is often hit first, so the search radius is set to 15arcmin
radius = '15' #arcmins. Vphas ccds are ~15 x 7.5 arcmin


#get the folder of g band CCD images 
#block_choice = raw_input('Block (a,b): ')



#create a dir to keep the panstarrs files in
pan_path = os.getcwd() + '/panstarrs'
if not os.path.exists(pan_path):
	os.makedirs(pan_path)
	print 'Created', pan_path
	downloaded_panstarrs = []
#if the directory does exists, see what files are in there
else:
	downloaded_panstarrs = glob.glob( pan_path+'/*')
		
	
	
for block_choice in ['a', 'b']:

	print 'Block ', block_choice

	#loop through the ccds and get the wcs coordinate of the central pixel
	ccds = glob.glob( os.getcwd() + '/vphas_*_ex/*g_*_'+block_choice+'/single/*_ccds/*.fit')

	for ccdnum in range(1,33):
		print 'CCD', ccdnum
	
		#path the downloaded data will be saved to. If it already exists in
		#the storage directory, continue
		savename = pan_path+'/block'+block_choice+'_ccd'+str(ccdnum)+'.csv'
		if savename in downloaded_panstarrs:
			print 'File already downloaded', savename
			print
			continue
	
		ccd_fpath = [fname for fname in ccds if '_ccd'+str(ccdnum)+'.fit' in fname]
	
		#open the ccd and get the wcs information
		ccd = fits.open( ccd_fpath[0] )
		hdr = ccd[1].header
		img = ccd[1].data
		ccd_wcs = wcs.WCS(hdr) 
	
		#central pixels as integars (skip the decimal point)
		cx = img.shape[0]/2
		cy = img.shape[1]/2
	
		#get the wcs of these central pixels
		ra , dec = ccd_wcs.wcs_pix2world(cx, cy, 1) 
		ccd.close()	
	

		# create URL with desired search parameters and
		#search the panstarrs query form
		url = "http://archive.stsci.edu/panstarrs/search.php?"  

		#add ra, dec and radius params
		url = url + "RA="+str(ra)+"&DEC="+str(dec)+"&radius="+radius

		#add, or rather remove, search limits
		#max_records: default = 2000. Limit may be set by available memory, but values > 25,000 may not 		complete.
		url = url + "&max_records=25000"
	
		#get coords back in degrees
		url = url + "&coordformat=dec"
	
		
	
		#output as csv
		url = url + "&outputformat=CSV"

		url = url + "&action=Search"
		print url

		# retrieve URL and  write results to filename
		urllib.urlretrieve(url,savename)
		print

print 'Panstarrs download complete'
print



#Combine the data from all the CCDS, removing repeats
print 'Combining data'

fits_cols = {}

for block_choice in ['a', 'b']:

	print 'Block', block_choice
	for ccdnum in range(1,33): 	
		print 'CCD', ccdnum
	
		savename = pan_path+'/block'+block_choice+'_ccd'+str(ccdnum)+'.csv'
		with open(savename, 'r') as f:
			reader = csv.reader(f)
			tab = [line for line in reader]
		

		colnames = tab[0]
		if colnames[0]=='no rows found':
			print 'No Panstarrs data available'
			continue
		
		types = tab[1]
		
		#get the column indicies for the number of Sloan band observations, and the magnitude columns
		ng_ind = [ ind for ind, colname in enumerate(colnames) if colname=='ng'][0]
		nr_ind = [ ind for ind, colname in enumerate(colnames) if colname=='nr'][0]
		ni_ind = [ ind for ind, colname in enumerate(colnames) if colname=='ni'][0]
		
		g_mag = [ ind for ind, colname in enumerate(colnames) if colname=='gMeanPSFMag'][0]
		r_mag = [ ind for ind, colname in enumerate(colnames) if colname=='rMeanPSFMag'][0]
		i_mag = [ ind for ind, colname in enumerate(colnames) if colname=='iMeanPSFMag'][0]
		


		tab = tab[2:]
			
		if len(tab)==0:
			print 'Table length of 0'
			continue
				
			
			
		#remove objects that haven't been observed at least once in g, r, i
		tab = [ line for line in tab if int(line[ng_ind])>0. and int(line[nr_ind])>0 and int(line[ni_ind])>0]
			
			
		#remove objects that don't have sensible magnitudes	
		tab = [ line for line in tab if float(line[g_mag])>0. and float(line[r_mag])>0 and float(line[i_mag])>0]
			
		
		#add to the dict
		for ind,name in enumerate(colnames):
			fits_cols[name] =  [line[ind] for line in tab]

	print
		
		

#create a astropy Table using the dictionary
all_cols = []
for ind, colname in enumerate(colnames):
	
	if types[ind]=='string':
		col = fits.Column(name=colname, format='20A', array = fits_cols[colname])
		
	else:
		#NEED MORE DECIMAL PLACES?
		col = fits.Column(name=colname, format='D', array = fits_cols[colname])
		
	all_cols.append(col)
	

hdu = fits.BinTableHDU.from_columns(all_cols)
hdu.writeto(fin_fpath)


#reopen the table to remove repeat entries
reopened = fits.open(fin_fpath, mode='update')
tab = reopened[1].data
new_tab = np.unique(tab)
reopened[1].data = new_tab
reopened.close()



print 'Final file saved: ', fin_fpath
print 

		































