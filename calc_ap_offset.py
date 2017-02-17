#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python


"""calculate offset between vphas apertures and apass (identical to first part of apass_colour_diagram.py)
Cross-match vphas aperture and apass using 12<vphas_mags<17. Use topcat to find all matches in a 5 arcsec radius and 
choose the brightest matched apass star as the true counterpart.
Create text file list of offsets for each aperture. These are equivalent to the apcor header values """

import os
from astropy.io import fits
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import glob
import math
import scipy.optimize
from datetime import datetime
import shutil

import make_lists
#change font size
import matplotlib
matplotlib.rcParams.update({'font.size': 25})



		

args = make_lists.get_vphas_num()
ex_path = os.getcwd()+'/vphas_'+args.vphas_num+'_ex'
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)

block_choice = raw_input('Choose block (a, b, c) : ')
if block_choice=='a': block = a_block
elif block_choice=='b': block = b_block
elif block_choice=='c': block = c_block


#apass file
#apass = os.getcwd() +'/apass_2.0.csv'
apass = os.getcwd() +'/apass_2.0_norepeats.csv'
#apass = os.getcwd() +'/apass_2.0_bestrepeats.csv'




#############################################################################################################
#remove any apass objects that have repeats
if apass == os.getcwd() +'/apass_2.0_norepeats.csv' and not os.path.exists(apass):
	print 'Removing any repeated objects in apass'
	apass_tab = np.recfromcsv(os.getcwd() +'/apass_2.0.csv')
	
	#new savepath
	new_apass = os.getcwd() +'/apass_2.0_norepeats.csv'	

	#remove entries that don't have g, r and i measurments
	titles = apass_tab.dtype.names
	#apass_tab = apass_tab[ apass_tab['sloan_g']!='NA' ]
	#apass_tab = apass_tab[ apass_tab['sloan_r']!='NA' ]
	#apass_tab = apass_tab[ apass_tab['sloan_i']!='NA' ]
	#print 'Removed incomplete entries'
	

	#remove entries with another entry within 0.0004 degrees
	newtab = []
	for i,line in enumerate(apass_tab):
		nearby = [j for j,row in enumerate(apass_tab) if i!=j and abs(float(row['radeg'])-float(line['radeg']))<0.0004 and abs(float(row['decdeg'])-float(line['decdeg']))<0.0004]
		
		if len(nearby)==0: #no repeats
			newtab.append(line)
	print 'Removed repeated entries'


	import csv
	with open(new_apass, 'wb') as f:
		_writer = csv.writer(f)
		_writer.writerow(titles)
		for line in newtab:
			_writer.writerow(line)
	print 'Written ', new_apass
	print



##############################################################################################################
#remove repeated objects and keep the "best" entry
if apass == os.getcwd() +'/apass_2.0_bestrepeats.csv' and not os.path.exists(apass):
	print 'Choosing best repeated entries in the apass catalogue'
	apass_tab = np.recfromcsv(os.getcwd() +'/apass_2.0.csv')
	new_apass = os.getcwd() +'/apass_2.0_bestrepeats.csv'	

	titles = apass_tab.dtype.names
	newtab = []
	for i,line in enumerate(apass_tab):
		print 'line', i
		nearby = [j for j,row in enumerate(apass_tab) if i!=j and abs(float(row['radeg'])-float(line['radeg']))<0.0004 and abs(float(row['decdeg'])-float(line['decdeg']))<0.0004]
	
		if len(nearby)==0:
			newtab.append(line)
			continue
		
		for j in nearby:
			if j<i: break #stops repeats eg. 3 matched with 5 then 5 matched with 3
			testline = apass_tab[j]
		
			#choose the one with the most observations
			if line['number_of_obs']>testline['number_of_obs']:
					newtab.append(line)
					break
			if testline['number_of_obs']>line['number_of_obs']:
					newtab.append(testline)
					break		
				
			#else choose the one with the fewest NAN measurments
			line_count = 0
			for val in line:
				if val=='NA': line_count+=1
				
			testline_count = 0
			for val in testline:
				if val=='NA': testline_count+=1
			
			if line_count>testline_count:
				newtab.append(testline)
				break	
				
			if testline_count>line_count:
				newtab.append(line)
				break
		
		
			#choose brightest V mag
			if line['johnson_v']!='NA' and testline['johnson_v']!='NA':
				if float(line['johnson_v'])<float(testline['johnson_v']): 
					newtab.append(line)
					break
				else:
					newtab.append(testline)
					break
				
			#else choose brightest B mag
			if line['johnson_b']!='NA' and testline['johnson_b']!='NA':
				if float(line['johnson_b'])<float(testline['johnson_b']): 
					newtab.append(line)
					break
				else:
					newtab.append(testline)
					break
				
			#else choose brightest g mag
			if line['sloan_g']!='NA' and testline['sloan_g']!='NA':
				if float(line['sloan_g'])<float(testline['sloan_g']): 
					newtab.append(line)
					break
				else:
					newtab.append(testline)
					break

				
			#else throw up error and keep neither
			#print 'Computer says no'
			#print line
			#print testline
			#print


	import csv
	with open(new_apass, 'wb') as f:
		_writer = csv.writer(f)	
		_writer.writerow(titles)
		for line in newtab:
			_writer.writerow(line)
	
	

###############################################################################################################


#Match vphas and apass objects.
#use topcat to match the apass to vphas ccd stars and save the resultant table
#Use magnitude using aperture 7 to remove bad stars

if block_choice=='a' or block_choice=='b':
	filternames = {'g':1, 'r_r':2, 'r_b':3, 'i':4}
elif block_choice=='c':
	filternames = {'g':1}
	
	
for filtername in filternames.keys():	
	catpath = glob.glob(block[filternames[filtername]]+'/catalogues/*cat.fits')
	print block[filternames[filtername]]
	catpath = catpath[0]
	
	#directory for apass-vphas merged catalogues
	fits_dirname = os.getcwd()+'/'+block_choice+'_'+filtername+'_fitsfiles'
	if not os.path.exists(fits_dirname):
		os.makedirs(fits_dirname)
	
	for ccdnum in range(1,33):	
	
		#merged apass+ccd filename
		apass_ccd_name = fits_dirname+'/apass_'+block_choice+'ccd'+str(ccdnum)+'.fits'
		#skip if the file already exists
		if os.path.exists(apass_ccd_name):
			continue
			
		print 'ccd ', str(ccdnum)
		
		#create a temporary fits file of the ccd catalogue to pipe into topcat
		temp_path = catpath[:-5]+'_ccd'+str(ccdnum)+'.fits'
	
		if os.path.exists(temp_path):
			os.remove(temp_path)
			
		cat = fits.open(catpath)
		all_hdus = []
		prihdr = cat[0].header
			
		all_hdus.append(fits.PrimaryHDU(header=prihdr))
		
		table = cat[ccdnum].data
		hdr = cat[ccdnum].header
		saturate = hdr['saturate']
		seeing = hdr['seeing'] #avergae fwhm in pixels
		cat.close()
	
		colnames = table.dtype.names
		colfields = table.dtype.fields
		
		#only use objects flagged as stars
		newtable = table[table['Classification']==-1] 
		new_cols = []
		for name in colnames:
			if 'Blank' in name: continue
			col = fits.Column(name=name, format=colfields[name][0], array=newtable[name])
			new_cols.append(col)
 
 		#add aperture 7 AB mags : use this to filter out bad stars
 		mags = [(-2.5*math.log10(line / hdr['exptime']) ) + hdr['nightzpt'] if line>0 else float('nan')  for line in table['Aper_flux_7'] ]

 		#remove vphas mags brighter than 12 and fainter than 20
 		mags = [line-float(hdr['APCOR7']) for line in mags] #apcor is the vphas catalogues aperture correction
 		mags = [line if 12.<line<17.0 else float('Nan') for line in mags ]
 		
 			
		mag_col = fits.Column(name='Aper_7_mag', format='D', array=mags)
		new_cols.append(mag_col)
		
		tbhdu = fits.BinTableHDU.from_columns(new_cols, header=hdr)
		all_hdus.append(tbhdu)
		tbhdulist = fits.HDUList(all_hdus)
		tbhdulist.writeto(temp_path, output_verify='silentfix')	
			

	

		#cross-match the apass and ccd file: use a large radius and keep all matches
		print
		print 'Merging apass and',block_choice,'block ccd', str(ccdnum)
		apass_match_radius = 5 #arcsec
		os.system('java -jar /mirror/scratch/hbarker/pkgs/jystilts.jar /home/hbarker/scripts/tmatch2.py {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(temp_path, apass, 'fits', True, 'csv', False, 'RA', 'DEC' , 'radeg', 'decdeg', apass_match_radius, '1and2', 'all', apass_ccd_name) )


		print
		print 'Removing repeat apass matches and keeping the brightest vphas star'
		matched = fits.open(apass_ccd_name, mode='update')
		table = matched[1].data
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
					
		for val in del_seq_nums:
			table = table[table['Sequence_number']!=val]


			
		
		
		
		"""We now have an apass - vphas ccd crossmatch list containing only the best matches. Now, 
		objects with poor photometry are removed."""
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




		"""cross-match the vphas ccd catalogue with the cross-matched file containing only the best apass-vphas
		photometry. Use this to calculate the distance between a matched object and any other nearby stars. If 
		there is a vphas star within 2x the seeing of the matched star, the matched star's vphas photometry will be
		contaminated, and so the match is removed."""

		#remove objects within the contamination distance
		print
		#temporary filename
		distances = fits_dirname+'/apass_'+block_choice+'ccd'+str(ccdnum)+'_distances.fits'
		os.system('java -jar /mirror/scratch/hbarker/pkgs/jystilts.jar /home/hbarker/scripts/tmatch2.py {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(apass_ccd_name, temp_path, 'fits', True, 'fits', True, 'RA', 'DEC' , 'RA', 'DEC', 15, '1and2', 'all', distances))
		
		distancefile = fits.open(distances)
		distancetab = distancefile[1].data
		distancefile.close()
		matched = fits.open(apass_ccd_name, mode='update')
		table = matched[1].data
		
		min_pixel_separation = 2 * seeing #seeing is fwhm in pixels
		min_separation = 0.218*min_pixel_separation #in arcsec, as topcat separation column is in arcsec
		
		#remove objects that are too close together by using the sequence numbers
		del_seq_nums = [line['Sequence_number_1'] for line in distancetab if 0.1<line['Separation']<min_separation and line['Sequence_number_1'] !=line['Sequence_number_2']  ]
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
corrections_dir = os.getcwd()+'/aperture_corrections'
if not os.path.exists(corrections_dir):
	os.makedirs(corrections_dir)	
	
#store details of any ccds with no vphas-apass matches	
no_matches = []	
		
print 'Calculating aperture corrections'
for filtername in filternames.keys():
	print filtername
	
	#directory containing the apass-vphas cross matched fits files
	fits_dirname = os.getcwd()+'/'+block_choice+'_'+filtername+'_fitsfiles'
	
	for apnum in range(2,8): #loop over apertures 2 to 7
		print 'Aperture', apnum
		
		#array to store aperture corrections. intersect_list = [ [ccdnum, value], ...] to keep track of cases with no aperture correction
		intersect_list = []
		corrections_path = corrections_dir+'/'+block_choice+'_'+filtername+'_aper'+str(apnum)+'_corrections.txt'
		if os.path.exists(corrections_path): 
			continue
			
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
			
			if len(table)==0:
				err_line = 'Block '+ str(block_choice) + ', filter: ' + str(filtername) + ', CCD: ' + str(ccdnum) + ', Ap: ' + str(apnum) + ', no matches '
				intersect_list.append( [ccdnum, float('nan') ] )
				error_list.append( [ccdnum, float('nan') ] )
				no_matches.append(err_line)
				continue
			
			
			#find vphas magntiudes NOTE: APASS and NIGHTZPT are in is AB
			vphas_AB = [(-2.5*math.log10(line / hdr['exptime']) ) + hdr['nightzpt'] if line>0 else float('nan')  for line in table['Aper_flux_'+str(apnum)] ]


			#apass mags
			if filtername=='r_r' or filtername=='r_b':
				apass_AB = table['Sloan_r']
			else:
				apass_AB = table['Sloan_'+filtername]		

				
			vphas_apass = [[line[0],line[1]] for line in zip(vphas_AB, apass_AB) if line[1]!='NA' and ~np.isnan(line[0])]	
			vphas_AB = [line[0] for line in vphas_apass]
			apass_AB = [float(line[1]) for line in vphas_apass]		
				
			#array: [difference between the apass and vphas magnitudes for the same object, vphas magnitude]
			difference_vphas = [[float(line[0])-line[1], float(line[0])] for line in zip(apass_AB, vphas_AB)]
		

			#plot graph including best fit line:
			x = [line[0] for line in difference_vphas]
			y= [line[1] for line in difference_vphas]
			
		
			#only use values in a well-behaved range ie. 13<g_vphas<16.  APASS saturates for mags <11 and becomes inaccurate >~17
			# have to swap the x and y axes
			newx = [line[1] for line in difference_vphas if 13<line[1]<16]
			newy = [line[0] for line in difference_vphas if 13<line[1]<16]
			
			allx = [line[1] for line in difference_vphas if 13<line[1]<16]
			ally = [line[0] for line in difference_vphas if 13<line[1]<16]
			
			
			#save an error if this filtering removes all the objects in the list
			if len(newx)==0 or len(newy)==0:
				print 'No suitable apass matches found for ccd', ccdnum
				err_line = 'Block '+ str(block_choice) + ', filter: ' + str(filtername) + ', CCD: ' + str(ccdnum) + 'Ap: ' + str(apnum) + ', no matches with suitable magnitudes '
				intersect_list.append( [ccdnum, float('nan') ] )
				error_list.append( [ccdnum, float('nan') ] )
				no_matches.append(err_line)
				continue
			

			#create function using numpy with a gradient of 1.0 for initial conditions
			coeffs = np.polyfit(newx, newy, 1)
			polynomial = np.poly1d(coeffs)
			init_intersect = polynomial[1]
			
			repeat = True
			while repeat==True:
			
				#use scipy to calculate the optimum y-intersection value
				errfunc = lambda c, y: c-y #ie. y = 0*x + c  -> y=c   ->  c-y=0
				intersection, success = scipy.optimize.leastsq(errfunc, init_intersect, args=newy)
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
				
				
			#save the intersection: intersection = g_apass - g_vphas  -> g_apass = g_vphas + intersection
			# in process_catalogues.py: g_vphas_true = g_vphas_measured - intersection
			
			intersect_list.append([ccdnum, 0-intersection]) #to keep consistent with how I'm adjusting mags	
			
			#error = sd of clipped mean	
			variances = [(val-intersection)**2 for val in newy]
			error = np.mean(variances)
			error_list.append([ccdnum, error])

				
			#final plot
			yrange = np.linspace(11,19)
			bestfit = [intersection for val in yrange]
			
			plt.plot(x,y ,'bo')
			plt.plot(newy, newx, 'ro')
			plt.plot(bestfit, yrange, 'r-')
			plt.ylim(11, 18.0)
			plt.gca().invert_yaxis() #so fainter stars at bottom of y axis
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
	



#delete working files if they're not needed 
 



#print out the list of any CCDs with no apass-vphas matches
if len(no_matches)>0:
	print
	print 'No APASS - VPHAS+ matches were found for: '
	for line in no_matches:
		print line
	print
	
	
	
	



