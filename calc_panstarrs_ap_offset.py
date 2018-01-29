#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python


"""calculate offset between vphas apertures and apass (identical to first part of apass_colour_diagram.py)
Cross-match vphas aperture and apass using 13<vphas_mags<16. Use topcat to find all matches in a 5 arcsec radius and 
choose the brightest matched apass star as the true counterpart.
Create text file list of offsets for each aperture. These are equivalent to the apcor header values """

import os
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import glob
import math
import scipy.optimize
import shutil
import argparse
import itertools
import sys
from astropy.table import Table, unique

import make_lists
#change font size
import matplotlib
matplotlib.rcParams.update({'font.size': 25})



#path to jystilts
jy_path = '/mirror/scratch/hbarker/pkgs/jystilts.jar'
#path to tmatch2.py. Is called as part of os, so must be the whole path
match2_fpath = '/home/hbarker/scripts/calibrate_vphas/tmatch2.py'





		
parser = argparse.ArgumentParser(description="Collect calibration inputs")
parser.add_argument('-v','--vphas_num', help="vphas pointing number", required=True)
parser.add_argument('-b','--block', help="vphas offset block", required=True)
args = parser.parse_args()


ex_path = os.getcwd()+'/vphas_'+args.vphas_num+'_ex'
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)







#block_choice = raw_input('Choose block (a, b, c) : ')
block_choice = args.block
if block_choice=='a': block = a_block
elif block_choice=='b': block = b_block
elif block_choice=='c': block = c_block

#panstarrs file
panstarrs = os.getcwd() +'/panstarrs.fits'
if not os.path.exists(panstarrs):
	print 'Could not find the reduced panstarrs file:', panstarrs
	sys.exit()
	

zptname = 'Nightzpt'	
	

#directory to store the lists of aperture corrections
corrections_dir = os.getcwd()+'/aperture_corrections'
if not os.path.exists(corrections_dir):
	os.makedirs(corrections_dir)	






#Match vphas and panstarrs objects.
#use topcat to match the panstarrs to vphas ccd stars and save the resultant table
#Use magnitude using aperture 7 to remove bad stars

if block_choice=='a' or block_choice=='b':
	filternames = {'g':1, 'r':2, 'r2':3, 'i':4} #u=0, NB=5

elif block_choice=='c':
	filternames = {'g':0} #NB=1
	
	
for filtername in filternames.keys():	
	print filtername
	
	
	#if the catalogue can't be found
	if block[filternames[filtername]] is None:
		print filtername, 'catalogue is None'
		sys.exit()
		

	catpath = glob.glob(block[filternames[filtername]] + '/catalogues/*cat_original.fits')
	if len(catpath)==0:
		catpath = glob.glob(block[filternames[filtername]] + '/catalogues/*cat.fits')
	
	if len(catpath)==0:
		print 'No catalogue'
		print block[filternames[filtername]]
		#sys.exit()
		raw_input('Press any key to continue')
		continue
		
		
	catpath = catpath[0]

	#directory for apass-vphas merged catalogues
	fits_dirname = os.getcwd()+'/'+block_choice+'_'+filtername+'_fitsfiles'
	print fits_dirname
	if not os.path.exists(fits_dirname):
		os.makedirs(fits_dirname)
	
	
	
	for ccdnum in range(1,33):	
	
	
		#merged apass+ccd filename
		panstarrs_ccd_name = fits_dirname+'/panstarrs_'+block_choice+'ccd'+str(ccdnum)+'.fits'
			
			
			
		#skip if the file already exists
		if os.path.exists(panstarrs_ccd_name):
			continue
			
			
		#filepath containing the aperture corrections
		#skip if it already exists	
		for apnum in range(2,5):
			corrections_path = corrections_dir+'/'+block_choice+'_'+filtername+'_aper'+str(apnum)+'_corrections.txt'
			if os.path.exists(corrections_path): 
				continue	
			
			
			
			
		print 'ccd ', str(ccdnum)
		
		
		
		
		#vphas catalogue files contain the catalogues for all 32 ccds
		#create a temporary fits file of the single ccd catalogue to pipe into topcat
		temp_path = catpath[:-5]+'_ccd'+str(ccdnum)+'.fits'
		print 'Creating temporary file: ', temp_path
	
		if os.path.exists(temp_path):
			os.remove(temp_path)
			
			
		cat = fits.open(catpath)
		all_hdus = []
		prihdr = cat[0].header
			
		all_hdus.append(fits.PrimaryHDU(header=prihdr))
		
		table = cat[ccdnum].data
		hdr = cat[ccdnum].header
		saturate = hdr['saturate'] #level at which the pixels saturate
		seeing = hdr['seeing'] #average fwhm in pixels
		cat.close()
	
		colnames = table.dtype.names
		colfields = table.dtype.fields
		
		#only use objects flagged as stars and remove 'Blank' columns
		newtable = table[table['Classification']==-1] 
		new_cols = []
		for name in colnames:
			if 'Blank' in name: continue
			if 'Aper_7_mag' in name: continue
			col = fits.Column(name=name, format=colfields[name][0], array=newtable[name])
			new_cols.append(col)
			
			
 
 		#add aperture 7 AB mags : use this to filter out bad stars
 		mags = [( -2.5*math.log10(line / hdr['exptime'] ) ) + hdr[zptname] if line>0. else float('nan')  for line in table['Aper_flux_7'] ]
 		
 		#apply the vphas catalogues aperture correction
 		#I find them consistent with the ones I calculate to ~0.001 mags
 		#mags = [line-float(hdr['APCOR7']) for line in mags] 

 		#remove vphas mags brighter than 13th and fainter than 21st to avoid poor vphas photometry
 		mags = [line if 13.<line<21.0 else float('Nan') for line in mags ]
 		
 			
		mag_col = fits.Column(name='Aper_7_mag', format='D', array=mags)
		new_cols.append(mag_col)
		
		
		tbhdu = fits.BinTableHDU.from_columns(new_cols, header=hdr)
		all_hdus.append(tbhdu)
		tbhdulist = fits.HDUList(all_hdus)
		tbhdulist.writeto(temp_path, output_verify='silentfix')	
		print 'Temporary vphas CCD file created'
		print	




	

		#cross-match the panstarrs and ccd file
		#panstarrs coords can be offset, so quite a large radius is needed
		print 'Merging panstarrs and',block_choice,'block ccd', str(ccdnum)
		
		#if filtername=='i': 
		#	match_radius = 1.0
		#else:
	#		match_radius = 3.0 #arcsec
			
		match_radius = 1.0 #arcsec

		command_line = 'java -jar '+jy_path+' '+match2_fpath+' {} {} {} {} {} {} {} {} {} {} {} {} {} {}'
		os.system(command_line.format(temp_path, panstarrs, 'fits', True, 'fits', False, 'RA', 'DEC' , 'raMean', 'decMean', match_radius, '1and2', 'all', panstarrs_ccd_name) )



		print
		print 'Removing repeat panstarrs matches'
		table = Table.read(panstarrs_ccd_name)
		print 'Table length:', len(table)
		
		#skip if there aren't any panstarrs-vphas matches
		if len(table)==0:
			os.remove(temp_path)
			print 'No matches: continuing'
			continue
			
		
		#remove any rows that aren't unique
		#newtab = unique(table, keys=['RA', 'DEC', 'raMean', 'decMean'])
		newtab = unique(table, keys='Sequence_number')
		
			
		
		print 'Number of unique rows:', len(newtab)
		
		#overwrite the previous file with the new one
		newtab.write(panstarrs_ccd_name, format='fits', overwrite=True)
		
		
		#clean up, removing temporary files
		os.remove(temp_path)

		


################################################################################################################
#plot panstarss vs vphas magntitudes for different apertures to calculate to the aperture to panstarrs correction

print
print
	
	
	
#store details of any ccds with no vphas-panstarrs matches	(because of poor coverage)
no_matches = []	



		
print 'Calculating aperture corrections'
for filtername in filternames.keys():


	#directory containing the panstarrs-vphas cross matched fits files
	fits_dirname = os.getcwd()+'/'+block_choice+'_'+filtername+'_fitsfiles'
	
	
	
	for apnum in range(2,8): #loop over apertures 2 to 7
		print 'Aperture', apnum
		
		
		#skip if the aperture corrections file already exists
		corrections_path = corrections_dir+'/'+block_choice+'_'+filtername+'_aper'+str(apnum)+'_corrections.txt'
		if os.path.exists(corrections_path): 
			continue
		
		
		#array to store aperture corrections
		#intersect_list = [ [ccdnum, value], ...] to keep track of cases with no aperture correction
		intersect_list = []
		

			
		#array for the error associated with the aperture correction
		error_list = []
		errors_path = corrections_dir+'/'+block_choice+'_'+filtername+'_aper'+str(apnum)+'_correction_errs.txt'



		for ccdnum in range(1,33):
		
			
			print 'ccd ', str(ccdnum)	
			
			#open the vphas catalogue to get header values
			catpath = glob.glob(block[filternames[filtername]]+'/catalogues/*cat.fits')
			print 'Catalogue:', catpath
			cat = fits.open(catpath[0])
			hdr = cat[ccdnum].header
			cat.close()


			#panstarrs - vphas crossmatched filename
			panstarrs_ccd_name = fits_dirname+'/panstarrs_'+block_choice+'ccd'+str(ccdnum)+'.fits'
			openfile = fits.open(panstarrs_ccd_name)
			table = openfile[1].data
			openfile.close()
			
			
			
			print 'Table length', len(table)
	
	
			
			#check to see if there are no panstarrs-vphas crossmatches
			if len(table)<5:
				err_line = 'Block '+ str(block_choice) + ', filter: ' + str(filtername) + ', CCD: ' + str(ccdnum) + ', Ap: ' + str(apnum) + ', no matches '
				intersect_list.append( [ccdnum, float('nan') ] )
				error_list.append( [ccdnum, float('nan') ] )
				no_matches.append(err_line)
				continue
			
			
			
			#calculate vphas magntiudes from the aperture count NOTE: APASS, APASSZPT and NIGHTZPT are in AB
			#Don't apply the aperture correction in the vphas header, we're calculating our own
			# If you do want to use it: mag = -2.5 * log10( counts * (1+apcor) / exptime ) + zpt
			vphas_AB = [ ( -2.5 * math.log10( line / hdr['exptime'] ) ) + hdr[zptname] if line>0 else float('nan')  for line in table['Aper_flux_'+str(apnum)] ]
			
			

			#Panstarrs seems to be in AB: https://arxiv.org/pdf/1203.0297.pdf
			#https://outerspace.stsci.edu/display/PANSTARRS/PS1+Database+object+and+detection+tables
			if filtername=='r2':
				pan_AB = table['rMeanPSFMag']
			else:
				pan_AB = table[filtername+'MeanPSFMag']		

		


			#remove any nan / NA entries in vphas/panstarrs
			vphas_pan = [[line[0],line[1]] for line in zip(vphas_AB, pan_AB) if line[1]!='NA' and ~np.isnan(line[0])]	
			vphas_AB = [line[0] for line in vphas_pan]
			pan_AB = [float(line[1]) for line in vphas_pan]	
			
			
				
			#array: [difference between the panstarrs and vphas magnitudes for the same object, panstarrs magnitude]
			difference_vphas = [[float(line[0])-line[1], float(line[0])] for line in zip(pan_AB, vphas_AB)]
		
		
			#remove any obvious mismatches; if the magnitude difference is more than 1.5 mags
			#difference_vphas = [ line for line in difference_vphas if abs(line[0])<1.5 ]
			
			
			
			#plot graph including best fit line:
			x = [line[0] for line in difference_vphas] # = panstarrs-vphas
			y= [line[1] for line in difference_vphas] # = panstarrs
			
					
		
			#only use values in a very well-behaved range ie. 13< mag <20
			# VPHAS is good up to ~21st mag, but the g band saturates around ~13th mag
			#Panstarrs saturation?
			
			
			# have to swap the x and y axes to fit the best fit line
			#newx, newy are used to fit the line of best fit
			if filtername=='g': low_limit = 13.0
			else: low_limit = 12.5
			
			high_limit = 20.0
			
			
			newx = [line[1] for line in difference_vphas if low_limit < line[1] < high_limit]  # = panstarrs
			newy = [line[0] for line in difference_vphas if low_limit < line[1] < high_limit] # = pan-vphas


			
			
						
			

			
			#flag an error if this filtering removes all the objects in the list
			if len(newx)==0 or len(newy)==0:
				print 'No suitable panstarrs matches found for ccd', ccdnum
				err_line = 'Block '+ str(block_choice) + ', filter: ' + str(filtername) + ', CCD: ' + str(ccdnum) + 'Ap: ' + str(apnum) + ', no matches with suitable magnitudes '
				intersect_list.append( [ccdnum, float('nan') ] )
				error_list.append( [ccdnum, float('nan') ] )
				no_matches.append(err_line)
				continue
			
			
			

			#create function using numpy with a gradient of 0 for initial conditions
			#init_intersect is an initial guess at the axis intercept
			coeffs = np.polyfit(newx, newy, 0)
			polynomial = np.poly1d(coeffs)
			init_intersect = polynomial[0]
		
			
			#test plot to make sure the initial conditions look ok
			#testy = [ init_intersect for val in newx]
			#plt.figure()
			#plt.plot(newx, newy, 'ro')
			#plt.plot( newx, testy, 'ko')
			#plt.show()

			
			intersection = None
			repeat = True
			while repeat==True:
			
				if intersection==None:
					intersection = init_intersect
			
				#use scipy to calculate the optimum y-intersection value for a fixed gradient
				errfunc = lambda c, y: c-y #ie. y = 0*x + c  -> y=c   ->  c-y=0
				intersection, success = scipy.optimize.leastsq(errfunc, intersection, args=newy)
				intersection = intersection[0] 
		
				#clip anything >3*median distance from best fit line and repeat the optimization if there are any points
				#further than this distance
				differences = [abs(line-intersection) for line in newy]
				limit = 3*np.median(differences)
				counter=0
		
				new_pairs = []
				for pair in zip(newx, newy):
					difference = abs(pair[1]-intersection)
					if difference>limit:
						counter+=1
					else:
						new_pairs.append(pair)
						
				newx = [line[0] for line in new_pairs]
				newy = [line[1] for line in new_pairs]
				
								
				if counter==0 : repeat=False
				


				
			#save the intersection: intersection = panstarrs - vphas  -> panstarrs = vphas + intersection
			# in process_catalogues.py: g_vphas_true = g_vphas_measured - intersection
			# so we save -1 *intersection
			intersect_list.append([ccdnum, 0-intersection]) 
			
			#calculate the error on the fit as the average distance to the line of best fit
			variance = [ (val-intersection)**2 for val in newy] 
			if len(variance)<=1:
				error = 'inf'
			else:
				st_dev =  math.sqrt( sum(variance) )
				error = st_dev / ( len(variance) - 1 )
			error_list.append([ccdnum, error])

				
			#final plot
			yrange = np.linspace(11,19)
			bestfit = [intersection for val in yrange]
			
			plt.plot(x,y ,'bo') #all points
			plt.plot(newy, newx, 'ro') #points used in the fit, remembering x and y have been swapped round
			plt.plot(bestfit, yrange, 'r-') #bestfit line
			plt.ylim(11, 18.0)
			plt.gca().invert_yaxis() #flip the y axis so fainter stars at bottom of y axis
			plt.subplots_adjust(left=None, bottom=0.13, right=None, top=None, wspace=None, hspace=None) #stop labels being clipped
		
					
			textline = 'x-intercept = ' + str(round(intersection,4))
			plt.annotate(textline, xy=(0.05, 0.90), xycoords='axes fraction', fontsize=24)
			
			if filtername=='r_r':
				graph_dirname = os.getcwd()+'/aperture_to_panstarrs/'+block_choice+'_r_aper'+str(apnum)+'_plots'
				catname1 = 'r_{Pan-STARRS}'
				catname2 = 'r_{VPHAS}'
			elif filtername=='r_b':
				graph_dirname = os.getcwd()+'/aperture_to_panstarrs/'+block_choice+'_r2_aper'+str(apnum)+'_plots'
				catname1 = '{r2}_{Pan-STARRS}'
				catname2 = '{r2}_{VPHAS}'
			else:
				graph_dirname = os.getcwd()+'/aperture_to_panstarrs/'+block_choice+'_'+str(filtername)+'_aper'+str(apnum)+'_plots'
				catname1 = str(filtername)+'_{Pan-STARRS}'
				catname2 = str(filtername)+'_{VPHAS}'
			
			#x = panstarrs - vphas
			plt.xlabel('$\mathregular{'+catname1+'}$ - $\mathregular{'+catname2+'}$', fontsize=32)
			plt.ylabel('$\mathregular{'+catname1+'}$', fontsize=32)	
			
			#plt.title('Aperture '+ str(apnum))

				
			if not os.path.exists(graph_dirname):
				os.makedirs(graph_dirname)
			savepath = graph_dirname + '/ccd'+str(ccdnum)+'.pdf'
			plt.savefig(savepath, bbox_inches='tight')
			#plt.show()
			plt.close()
			
	
	
		#write text file with aperture corrections
		txtfile = open(corrections_path, 'w')
		for line in intersect_list:
			#txtline = ccdnum   aperture correction
			txtline = str(line[0]) + '     ' + str(line[1]) + '\n'	
			txtfile.write(txtline)
		txtfile.close()
	


		#write file with aperture correction errors
		txtfile = open(errors_path, 'w')
		for line in error_list:
			line = str(line[0]) + '     ' + str(line[1]) + '\n'	
			txtfile.write(line)
		txtfile.close()
	





#print out the list of any CCDs with no panstarrs-vphas matches
if len(no_matches)>0:
	print
	print 'No Pan-STARRS - VPHAS+ matches were found for: '
	for line in no_matches:
		print line
	print
	
	
	
	



