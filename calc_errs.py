#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""calculate errros from photon count, aperture correction and zeropoint errors and add them to the raw catalogue (ie. the same one combine_cats.py/process_cats.py edits"""

import glob
from astropy.io import fits
import os
import numpy as np
import math

import make_lists

def append_table(table, name, arr, d_type):
    arr = np.asarray(arr)
    dtype = d_type
    newdtype = np.dtype(table.dtype.descr + [(name, d_type)])
    newtable = np.empty(table.shape, dtype=newdtype)
    for field in table.dtype.fields:
        newtable[field] = table[field]
    newtable[name] = arr
    return newtable
    
 
#vega to AB magnitude conversions
u_conv = 2.5*math.log((3631.0/1550.81) ,10)
g_conv = 2.5*math.log((3631.0/3960.53) ,10)
r_conv = 2.5*math.log((3631.0/3094.68) ,10)
i_conv = 2.5*math.log((3631.0/2563.84) ,10)
NBa_conv = 2.5*math.log((3631.0/2903.53) ,10)
NBb_conv = 2.5*math.log((3631.0/2931.14) ,10)
NBc_conv = 2.5*math.log((3631.0/2929.13) ,10)
NBd_conv = 2.5*math.log((3631.0/2669.32) ,10)


args = make_lists.get_vphas_num()
ex_path = os.getcwd()+'/vphas_'+args.vphas_num+'_ex'
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)

block_choice = raw_input('Choose block (a, b, c) : ')
if block_choice=='a': block = a_block
elif block_choice=='b': block = b_block
elif block_choice=='c': block = c_block

if block_choice=='a' or block_choice=='b':
	filternames = {'u':0, 'g':1,'r_r':2,'r_b':3,'i':4, 'NB':5}
elif block_choice=='c':
	filternames = {'g':1}

for i,filtername in enumerate(filternames):


	print filtername
	catpath = glob.glob(block[filternames[filtername]]+'/catalogues/*_fix_cat.fits')
	opencat = fits.open(catpath[0], mode='update')
	for ccdnum in range(1,33):

	
		print'ccd', ccdnum
		hdr = opencat[ccdnum].header
		table = opencat[ccdnum].data
		
		
		for apnum in range(2, 8):

		
			#skip anything but aperture 2,3,4,5 for u
			if filtername=='u' and apnum>5:
					continue
			
			
			ap_name = 'Aper_flux_'+str(apnum)
			print ap_name
			if ap_name+'_corr' not in table.dtype.names:
				continue
			mags = table[ap_name+'_corr']
			
			
			
			#exposure_time
			exp_t = table['Exp_time']
			
			#counts and error
			counts = table[ap_name]
			counts_err = table[ap_name+'_err']
			
			
			if filtername=='NB':
				highest_counts = np.add(counts, counts_err)
				low_mag = [ -2.5*math.log10( line[0] /line[1] ) + hdr['nightzpt'] if line[0]>0 else float('nan') for line in zip(highest_counts, exp_t) ]
				
				lowest_counts = np.subtract(counts, counts_err)
				high_mag = [ -2.5*math.log10( line[0] /line[1] ) + hdr['nightzpt'] if line[0]>0 else float('nan') for line in zip(lowest_counts, exp_t) ]
				
				#convert to vega:
	                	if ccdnum in range(0,9): #A
	                	        conv = NBa_conv
	                	elif ccdnum in range(9,17): #D
	                        	conv = NBd_conv
	                	elif ccdnum in range(17, 25):#B
	                        	conv = NBb_conv
	                	elif ccdnum in range(25, 33): #C
	                        	conv = NBc_conv
	                        	
	                       	low_mag = [line-conv for line in low_mag]
	              		high_mag = [line-conv for line in high_mag]
	              		
	              		
	              		#append column to table
				if not ap_name+'_upper_lim' in table.dtype.names:
					table = append_table(table, ap_name+'_upper_lim', high_mag, '>f4')
				else:
					table[ap_name+'_upper_lim'] = high_mag
			
				if not ap_name+'_lower_lim' in table.dtype.names:
					table = append_table(table, ap_name+'_lower_lim', low_mag, '>f4')
				else:
					table[ap_name+'_lower_lim'] = low_mag
	
				continue
	                        
	                        	
			
			
			#aperture correction and error in magnitudes
			elif filtername!='u':
				apcor_path = os.getcwd()+'/aperture_corrections/'+block_choice+'_'+filtername+'_aper'+str(apnum)+'_corrections.txt'
				with open(apcor_path, 'r') as f:
					for line in f:
						line = line.strip().split()
						if int(line[0])==ccdnum:
							apcor = float(line[1])

				apcor_err_path = os.getcwd()+'/aperture_corrections/'+block_choice+'_'+filtername+'_aper'+str(apnum)+'_correction_errs.txt'
				with open(apcor_err_path, 'r') as f:
					for line in f:
						line = line.strip().split()
						if int(line[0])==ccdnum:
							ap_error = float(line[1])
							
			#else if filtername==u		
			else:
				apcor_path = os.getcwd()+'/u_corrections_'+block_choice+'.txt'
				if not os.path.exists(apcor_path):
					print 'u corrections could not be found'
					print apcor_path
					
				with open(apcor_path, 'r') as f:
					for line in f:
						line = line.strip().split()
						if line[0]==str(apnum):
							apcor = float(line[1])
							ap_error = float(line[2])

			
			
			
			#convert aperture correction to counts	
			corrected_counts = [line[0]*10**((hdr['nightzpt']-line[1])/2.5) for line in zip(exp_t, mags)]	
			
			apcor_counts = np.subtract(corrected_counts, counts)


			#lowest_mag = brightest ie. most counts
			low_mag = [ ( -2.5*math.log10( (line[0]+line[1]) /line[2] ) + hdr['nightzpt'] - apcor - ap_error ) if line[0]+line[1]>0 else float('nan') for line in zip(counts, counts_err, exp_t) ]
			highest_counts = [line[0]*10**((hdr['nightzpt']-line[1])/2.5) for line in zip(exp_t, low_mag)]
			ap_upper_err_counts = np.subtract(highest_counts, corrected_counts)

			
			#highest mag = fewest counts
			high_mag = [ ( -2.5*math.log10( (line[0]-line[1]) /line[2] ) + hdr['nightzpt'] - apcor + ap_error ) if line[0]-line[1]>0 else float('nan') for line in zip(counts, counts_err, exp_t) ]
			lowest_counts = [line[0]*10**((hdr['nightzpt']-line[1])/2.5) for line in zip(exp_t, high_mag)]
			ap_lower_err_counts = np.subtract(corrected_counts, lowest_counts)


			#these upper and lower limits are in AB magnitudes. Need to convert to vega
			if filtername == 'r_r' or filtername == 'r_b': conv = r_conv
	            	if filtername == 'i': conv = i_conv
	                if filtername == 'u': conv = u_conv
	                if filtername == 'g': conv = g_conv
	              	low_mag = [line-conv for line in low_mag]
	              	high_mag = [line-conv for line in high_mag]
			
			""" #For checking mags
			mag = [ ( -2.5*math.log10( line[0] /line[1] ) + hdr['nightzpt'] -apcor ) if line[0]+line[1]>0 else float('nan') for line in zip(counts, exp_t) ]
			mag = [line-conv for line in mag]
			
			print mag
			print low_mag
			print high_mag
			print
			print mag[0]-low_mag[0]
			print mag[0]-high_mag[0]
			avg = math.sqrt( (mag[0]-low_mag[0])**2 + mag[0]-high_mag[0]**2)
			raw_input('')
			"""
	
	
			#append column to table
			if not ap_name+'_upper_lim' in table.dtype.names:
				table = append_table(table, ap_name+'_upper_lim', high_mag, '>f4')
			else:
				table[ap_name+'_upper_lim'] = high_mag
			
			if not ap_name+'_lower_lim' in table.dtype.names:
				table = append_table(table, ap_name+'_lower_lim', low_mag, '>f4')
			else:
				table[ap_name+'_lower_lim'] = low_mag
			
			
		opencat[ccdnum].data = table
	print 'Catalogue updated with error values'
	print
	opencat.close()
			
			
			
		

	
	
