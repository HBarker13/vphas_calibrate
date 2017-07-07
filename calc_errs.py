#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""calculate errros from photon count, aperture correction and zeropoint errors and add them to the raw catalogue (ie. the same one combine_cats.py/process_cats.py edits"""

import glob
from astropy.io import fits
import os
import numpy as np
import math
import argparse

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



if block_choice=='a' or block_choice=='b':
	filternames = {'u':0, 'g':1,'r':2,'r2':3,'i':4, 'NB':5}
elif block_choice=='c':
	filternames = {'g':0}
	


for i,filtername in enumerate(filternames):
	
	
	print filtername
	catpath = glob.glob( block[filternames[filtername]] + '/catalogues/*_fix_cat.fits')
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
				
				
			
			#Not sure this is correct. If you're using Halpha mags, I wouldn't trust this
			if filtername=='NB':
				continue
				"""
				highest_counts = np.add(counts, counts_err)
				low_mag = [ -2.5*math.log10( line[0] /line[1] ) + hdr['apasszpt'] if line[0]>0 else float('nan') for line in zip(highest_counts, table['Exp_time']) ]
				
				lowest_counts = np.subtract(counts, counts_err)
				high_mag = [ -2.5*math.log10( line[0] /line[1] ) + hdr['apasszpt'] if line[0]>0 else float('nan') for line in zip(lowest_counts, table['Exp_time']) ]
				
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
				"""
			
	                    
	                        	
			
			
			#read in aperture corrections and errors in magnitudes
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
							
							
		
	



			
			
			#mag = -2.5* log10( counts/ exp_t ) + zpt
			#counts =  exp_t * 10 ** ((zpt - mag)/2.5 )
			
			#calculate the zero point of the corrected vega magnitudes as Nightzpt (in AB) - vega to AB conversion
			#It doesn't matter that these aren't strictly correct for the u band, the aperture correction will adjust for this.
			if filtername=='u':
				zpt = [ line['Apasszpt'] - u_conv for line in table]
			if filtername=='g':
				zpt = [ line['Apasszpt'] - g_conv for line in table]
			if filtername=='r':
				zpt = [ line['Apasszpt'] - r_conv for line in table]
			if filtername=='r2':
				zpt = [ line['Apasszpt'] - r_conv for line in table]
			if filtername=='i':
				zpt = [ line['Apasszpt'] - i_conv for line in table]
			
			
			
			
			
			#combine the photon count error and aperture correction error (in counts) in quadrature
			
			#raw counts provided by the vphas catalogue
			counts = table[ ap_name ]
			counts_err = table[ap_name+'_err']

			#vphas header data used to calculate mags
			exp_time = table['exp_time']
			zpt = table['Apasszpt']
			
	
			#counts corresponding to the corrected magnitude
			corr_mag = table[ ap_name + '_corr']
			corr_counts = [ t*10**( (line - z) /-2.5) for t, z, line in zip(exp_time, zpt, corr_mag) ]
			
			
			
			
			#upper_mag = corrected_mag + aperture correction error = fewest counts
			upper_mag = [ line + ap_error for line in corr_mag ]
			upper_mag_counts = [ t*10**( (line - z )/-2.5) for t, z, line in zip(exp_time, zpt, upper_mag) ]
					
			#the number of counts corresponding to the aperture correction error			
			upper_ap_err_counts = np.subtract( corr_counts, upper_mag_counts )
			
			
			#combine the error from the counts, and the aperture correction error
			upper_tot_count_err = [ math.sqrt( line[0]**2 + line[1]**2 ) for line in zip(counts_err, upper_ap_err_counts) ]
			
			#upper mag = largest number = fewest counts
			upper_mag_lim_counts = np.subtract(corr_counts, upper_tot_count_err )
			high_mag = [ -2.5*math.log10( line / t ) + z if line>0 else float('nan') for t, z, line in zip(exp_time, zpt, upper_mag_lim_counts) ]
			
			
			
			
			
			#lower_mag = corrected_mag - aperture correction error
			lower_mag = [ line - ap_error for line in corr_mag ]
			lower_mag_counts = [ t*10**( (line - z )/-2.5) for t, z, line in zip(exp_time, zpt, lower_mag) ]
			
			#the number of counts corresponding to the aperture correction error			
			lower_ap_err_counts = np.subtract( lower_mag_counts,  corr_counts )
			
			#combine the error from the counts, and the aperture correction error
			lower_tot_count_err = [ math.sqrt( line[0]**2 + line[1]**2 ) for line in zip(counts_err, lower_ap_err_counts) ]
			
			#lower mag = smallest number = most counts
			lower_mag_lim_counts = np.add(corr_counts, lower_tot_count_err )
			low_mag = [ -2.5*math.log10( line / t ) + z if line>0 else float('nan') for t, z, line in zip(exp_time, zpt, lower_mag_lim_counts) ]
			
			

			



	
			#append column to table
			#upper_lim = highest magnitude number = fewest counts
			if ap_name+'_upper_lim' not in table.dtype.names:
				table = append_table(table, ap_name+'_upper_lim', high_mag, '>f4')
			else:
				table[ap_name+'_upper_lim'] = high_mag
			
			
			#lower_lim = lowest magnitude limit = highest counts
			if ap_name+'_lower_lim' not in table.dtype.names:
				table = append_table(table, ap_name+'_lower_lim', low_mag, '>f4')
			else:
				table[ap_name+'_lower_lim'] = low_mag
			
			opencat[ccdnum].data = table
	print 'Catalogue updated with error values'
	print
	opencat.close()
			
			
			
		

	
	
