#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python


#Read in the apass file in the directory. Remove /merge repeat measurments with different coordiantes. On;y keep entries with 2 or more observations



import os
import math
import numpy as np
import csv
import glob
import sys



#name of the new file to save to
savepath = os.getcwd()+'/apass_reduced.csv'


if os.path.exists(savepath):
	print 'File already exists!'
	print savepath
	sys.exit()




#read in the apass file (downlaoded from the apass webpage in csv format)
apass = [ fpath for fpath in glob.glob( os.getcwd() +'/*apass*.csv') if 'reduced' not in fpath]

if len(apass)==0:
	print 'No apass files found'
	print os.getcwd()
	raw_input('Paused')
	 




else:

	apass = apass[0]

	apass_tab = np.recfromcsv( apass )
	colnames = apass_tab.dtype.names

	#remove entries that don't have g, r and i measurments
	apass_tab = apass_tab[ apass_tab['sloan_g']!='NA' ]
	apass_tab = apass_tab[ apass_tab['sloan_r']!='NA' ]
	apass_tab = apass_tab[ apass_tab['sloan_i']!='NA' ]
	print 'Removed incomplete entries'
	

	#change dtype from stings to float. 'NA' is a string.
	new_dtype = [ (name, np.float64) if name!='johnson_v' and name!='verr' and name!='johnson_b' and name!='b_err' else (name, '|S6') for name in colnames ]
	apass_tab = apass_tab.astype( new_dtype )

	
	#remove objects with negative errors - no idea why they're in there
	apass_tab = apass_tab[ apass_tab['gerr']>0.0 ]
	apass_tab = apass_tab[ apass_tab['r_err']>0.0 ]
	apass_tab = apass_tab[ apass_tab['i_err']>0.0 ]
	print 'Removed entries with negative errors'
	
	
	
	#remove objects with errors < 0.1 mags
	apass_tab = apass_tab[ apass_tab['gerr']<1.0 ]
	apass_tab = apass_tab[ apass_tab['r_err']<1.0 ]
	apass_tab = apass_tab[ apass_tab['i_err']<1.0 ]
	print 'Removed entries with errors > 0.1 magnitudes'
	
	
	
	
	newtab = []
	#I've found apass photometry to be good within ~1 -1.5 arcseconds
	for line in apass_tab:
		ra = line['radeg']
		dec = line['decdeg']
		
		#want to match objects within 1.5 arcsecond radius
		#so within sqrt(( difference**2) / 2 ) along ra and dec
		difference = 1.5/3600.
		difference = math.sqrt( difference**2 / 2 )

		
		#matches within the differene distance
		matches = apass_tab[ abs(apass_tab['radeg']-ra) < difference ]
		matches = matches[ abs(matches['decdeg']-dec) < difference ]

		if len(matches)>1:
		
			#if the g, r and i magnitudes match to within 0.1  magnitudes, it is likely the same star
			mag_diff = 0.1
			newline = matches[0]
			ref_g = newline['sloan_g']
			ref_r = newline['sloan_r']
			ref_i = newline['sloan_i']

			for i in range(1, len(matches)): #skip comparing matches[0] to itself

				if abs(matches[i]['sloan_g']-ref_g)<mag_diff and abs(matches[i]['sloan_r']-ref_r)<mag_diff and abs(matches[i]['sloan_i']-ref_i)<mag_diff:
				
					
					newline['number_of_obs'] = len(matches)
			
					#the magnitude difference is small enough its probably not worth converting to counts
					g_avg = np.mean(matches['sloan_g'])
					r_avg = np.mean(matches['sloan_r'])
					i_avg = np.mean(matches['sloan_i'])
			
			
					#I don't use the errors so I won't bother updating the values
					newtab.append(newline)
					continue
		
		elif len(matches)==1:  #matched to self
			newtab.append(line)
			
	
	#remove objects with only one observation, even after combining those with erronous coordinates
	newtab = [line for line in newtab if line[4]>1.0 ] #line[4]=number_of_obs		


	#write the new apass data to the savepath
	with open(savepath, 'wb') as f:
		_writer = csv.writer(f)
		_writer.writerow(colnames)
		for line in newtab:
			_writer.writerow(line)
	print 'Written ', savepath
	print


	
	



