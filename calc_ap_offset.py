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

# reduced apass file
apass = os.getcwd() +'/apass_reduced.csv'
if not os.path.exists(apass):
	print 'Could not find the reduced apass file:', apass
	print 'Try reduce_apass.py'
	sys.exit()
	

#directory to store the lists of aperture corrections
corrections_dir = os.getcwd()+'/aperture_corrections'
if not os.path.exists(corrections_dir):
	os.makedirs(corrections_dir)	






#Match vphas and apass objects.
#use topcat to match the apass to vphas ccd stars and save the resultant table
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
		apass_ccd_name = fits_dirname+'/apass_'+block_choice+'ccd'+str(ccdnum)+'.fits'
			
			
			
		#skip if the file already exists
		if os.path.exists(apass_ccd_name):
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
			col = fits.Column(name=name, format=colfields[name][0], array=newtable[name])
			new_cols.append(col)
			
			
 
 		#add aperture 7 AB mags : use this to filter out bad stars
 		mags = [( -2.5*math.log10(line / hdr['exptime'] ) ) + hdr['apasszpt'] if line>0. else float('nan')  for line in table['Aper_flux_7'] ]
 		
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
			




	

		#cross-match the apass and ccd file: use a large radius and keep all matches
		#need to call the python script with stilts imported - doing it here breaks other things
		print 'Merging apass and',block_choice,'block ccd', str(ccdnum)
		apass_match_radius = 5 #arcsec

		command_line = 'java -jar '+jy_path+' '+match2_fpath+' {} {} {} {} {} {} {} {} {} {} {} {} {} {}'
		os.system(command_line.format(temp_path, apass, 'fits', True, 'csv', False, 'RA', 'DEC' , 'radeg', 'decdeg', apass_match_radius, '1and2', 'all', apass_ccd_name) )


		print
		print 'Removing repeat apass matches and keeping the brightest vphas star'
		matched = fits.open(apass_ccd_name, mode='update')
		table = matched[1].data
		
	
		
		#sequence numbers are unique within a ccd, so use them to identify objects that need removing
		del_seq_nums = set()
		for line in table:
			#find multiple matches to same apass object
			temptable = table[table['radeg']==line['radeg']]
			
			if len(temptable)==1: 
				continue #only matched to self

			#keep the brightest vphas object as the counterpart
			maxmag = min(temptable['Aper_7_mag']) #np.min includes nan, min doesn't
			for repeat in temptable:
				if repeat['Aper_7_mag']!=maxmag:
					del_seq_nums.add(repeat['Sequence_number'])
					
			
			
		#remove apass-vphas matches that are more than 1.5 arcseconds apart
		#This shouldn't be needed as we already match to the brightest star in vphas 
		#ie. If there was a brighter star in the vicinity of the apass star, that would be the apass star	
		for line in table:
			if line['Separation']>2.5:
				del_seq_nums.add( line['Sequence_number'])			
					
					
		for val in del_seq_nums:
			table = table[table['Sequence_number']!=val]


		
		
		
		
		#We have an apass - vphas ccd crossmatch list containing only the best matches. Now, 
		#objects with poor photometry are removed.
		#remove matches if the vphas object has an error bit flag
		table = table[table['error_bit_flag']==0]
		
		#remove vphas saturated objects
		del_seq_nums = set()
		for line in table:
			if line['peak_height']+line['sky_level'] > saturate: #saturate value from header
				del_seq_nums.add(line['Sequence_number'])
		for val in del_seq_nums:
			table = table[table['Sequence_number']!=val]
				
		matched[1].data = table
		matched.close()		
		




		#cross-match the vphas ccd catalogue with the cross-matched file containing the best apass-vphas
		#photometry. Use this to calculate the distance between a matched object and any other nearby stars. If 
		#there is a vphas star within 2x the seeing of the matched star, the matched star's vphas photometry will be
		#contaminated, and so the match is removed
		

		#temporary filename
		distances = fits_dirname+'/apass_'+block_choice+'ccd'+str(ccdnum)+'_distances.fits'
		
		command_line = 'java -jar '+jy_path+' '+match2_fpath+' {} {} {} {} {} {} {} {} {} {} {} {} {} {}'
		os.system(command_line.format(apass_ccd_name, temp_path, 'fits', True, 'fits', True, 'RA', 'DEC' , 'RA', 'DEC', 15, '1and2', 'all', distances))
		

		distancefile = fits.open(distances)
		distancetab = distancefile[1].data
		distancefile.close()
		matched = fits.open(apass_ccd_name, mode='update')
		table = matched[1].data
		
		
		min_pixel_separation = 2 * seeing #seeing is fwhm in pixels
		#OmegaCAM has 0.218 arcsec per pixel
		min_separation = 0.218 * min_pixel_separation #in arcsec, as topcat separation column is in arcsec
		

		
		
		#remove objects that are too close together by using the sequence numbers
		#line['Sequence_number_1'] != line['Sequence_number_2'] means objects aren't matched to themselves
		del_seq_nums = [line['Sequence_number_1'] for line in distancetab if line['Separation']<min_separation and line['Sequence_number_1'] != line['Sequence_number_2']  ]
		
		
		for val in del_seq_nums:
			table = table[table['Sequence_number']!=val]
		
		
		matched[1].data = table
		matched.close()			


		#clean up, removing temporary files
		os.remove(temp_path)
		os.remove(distances)
		


################################################################################################################
#plot apass vs vphas magntitudes for different apertures to calculate to the aperture to apass correction

print
	
	
	
#store details of any ccds with no vphas-apass matches	(because of poor apass coverage)
no_matches = []	



		
print 'Calculating aperture corrections'
for filtername in filternames.keys():
	print filtername


	#directory containing the apass-vphas cross matched fits files
	fits_dirname = os.getcwd()+'/'+block_choice+'_'+filtername+'_fitsfiles'
	
	
	
	for apnum in range(2,8): #loop over apertures 2 to 7
		print 'Aperture', apnum
		
		
		#skip if the aperture corrections file already exists
		corrections_path = corrections_dir+'/'+block_choice+'_'+filtername+'_aper'+str(apnum)+'_corrections.txt'
		#if os.path.exists(corrections_path): 
		#	continue
		
		
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
			cat = fits.open(catpath[0])
			hdr = cat[ccdnum].header
			cat.close()


			#apass - vphas crossmatched filename
			apass_ccd_name = fits_dirname+'/apass_'+block_choice+'ccd'+str(ccdnum)+'.fits'
			openfile = fits.open(apass_ccd_name)
			table = openfile[1].data
			openfile.close()
	
	
			
			#check to see if there are no apass-vphas crossmatches
			if len(table)<5:
				err_line = 'Block '+ str(block_choice) + ', filter: ' + str(filtername) + ', CCD: ' + str(ccdnum) + ', Ap: ' + str(apnum) + ', no matches '
				intersect_list.append( [ccdnum, float('nan') ] )
				error_list.append( [ccdnum, float('nan') ] )
				no_matches.append(err_line)
				continue
			
			
			
			#calculate vphas magntiudes from the aperture count NOTE: APASS and APASSZPT are in is AB
			#Don't apply the aperture correction in the vphas header, we're calculating our own
			# If you do want to use it: mag = -2.5 * log10( counts * (1+apcor) / exptime ) + zpt
			vphas_AB = [ ( -2.5 * math.log10( line / hdr['exptime'] ) ) + hdr['apasszpt'] if line>0 else float('nan')  for line in table['Aper_flux_'+str(apnum)] ]



			#apass mags
			if filtername=='r2':
				apass_AB = table['Sloan_r']
			else:
				apass_AB = table['Sloan_'+filtername]		



			#remove any nan / NA entries in vphas/apass
			vphas_apass = [[line[0],line[1]] for line in zip(vphas_AB, apass_AB) if line[1]!='NA' and ~np.isnan(line[0])]	
			vphas_AB = [line[0] for line in vphas_apass]
			apass_AB = [float(line[1]) for line in vphas_apass]		
				
				
			#array: [difference between the apass and vphas magnitudes for the same object, apass magnitude]
			difference_vphas = [[float(line[0])-line[1], float(line[0])] for line in zip(apass_AB, vphas_AB)]
		
		
			#remove any obvious mismatches; if the magnitude difference is more than 1.5 mags
			difference_vphas = [ line for line in difference_vphas if abs(line[0])<1.5 ]		




			
			#plot graph including best fit line:
			x = [line[0] for line in difference_vphas] # = apass-vphas
			y= [line[1] for line in difference_vphas] # = apass
			
					
		
			#only use values in a very well-behaved range ie. 13< mag <16. 
			# APASS saturates for mags <11 and becomes inaccurate >~16
			# VPHAS is good up to ~21st mag, but the g band saturates around ~13th mag
			
			
			# have to swap the x and y axes to fit the best fit line
			#newx, newy are used to fit the line of best fit
			if filtername=='g': low_limit = 13.0
			else: low_limit = 12.5
			
			high_limit = 17.0
			
			
			newx = [line[1] for line in difference_vphas if low_limit < line[1] < high_limit]  # = apass
			newy = [line[0] for line in difference_vphas if low_limit < line[1] < high_limit] # = apass-vphas


			
			
						
			

			
			#flag an error if this filtering removes all the objects in the list
			if len(newx)==0 or len(newy)==0:
				print 'No suitable apass matches found for ccd', ccdnum
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
				


				
			#save the intersection: intersection = apass - vphas  -> apass = vphas + intersection
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
				graph_dirname = os.getcwd()+'/aperture_to_apass/'+block_choice+'_r_aper'+str(apnum)+'_plots'
				catname1 = 'r_{APASS}'
				catname2 = 'r_{VPHAS}'
			elif filtername=='r_b':
				graph_dirname = os.getcwd()+'/aperture_to_apass/'+block_choice+'_r2_aper'+str(apnum)+'_plots'
				catname1 = '{r2}_{APASS}'
				catname2 = '{r2}_{VPHAS}'
			else:
				graph_dirname = os.getcwd()+'/aperture_to_apass/'+block_choice+'_'+str(filtername)+'_aper'+str(apnum)+'_plots'
				catname1 = str(filtername)+'_{APASS}'
				catname2 = str(filtername)+'_{VPHAS}'
			
			#x = apass - vphas
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
	





#print out the list of any CCDs with no apass-vphas matches
if len(no_matches)>0:
	print
	print 'No APASS - VPHAS+ matches were found for: '
	for line in no_matches:
		print line
	print
	
	
	
	



