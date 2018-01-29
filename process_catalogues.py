#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#Find catalogue fits files in expath and create a fits file containing the entries of all 32 ccds
#Adds some useful columns: filter, directory+ccd, exposure_time, zpt, mags
#New files are saved in directory in expath


import os
import glob
import subprocess as sp
from astropy.io import fits
from datetime import datetime
import numpy as np
from astropy.table import Column
import warnings
import math
import shutil

import make_lists



args = make_lists.get_vphas_num()
ex_path = os.getcwd() +  '/vphas_'+args.vphas_num+'_ex'
zptname = 'Nightzpt'
start = datetime.now()



#vega to AB magnitude conversions
u_conv = 2.5*math.log((3631.0/1550.81) ,10)
g_conv = 2.5*math.log((3631.0/3960.53) ,10)
r_conv = 2.5*math.log((3631.0/3094.68) ,10)
i_conv = 2.5*math.log((3631.0/2563.84) ,10)
NBa_conv = 2.5*math.log((3631.0/2903.53) ,10)
NBb_conv = 2.5*math.log((3631.0/2931.14) ,10)
NBc_conv = 2.5*math.log((3631.0/2929.13) ,10)
NBd_conv = 2.5*math.log((3631.0/2669.32) ,10)



#Loop though a, b, and c block in order to make everything easier to follow
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)

a_catnames = [] #in order u,g, r, r2, i,NB
for dirname in a_block:
	catname = glob.glob(dirname + '/catalogues/*cat.fits')
	if len(catname)!=0:
		a_catnames.append(catname[0])
		
		b_catnames = []
for dirname in b_block:
	catname = glob.glob(dirname + '/catalogues/*cat.fits')
	if len(catname)!=0:
		b_catnames.append(catname[0])
			
if None not in c_block:	
	c_catnames = []
	for dirname in c_block:
		catname = glob.glob(dirname + '/catalogues/*cat.fits')
		if len(catname)!=0:
			c_catnames.append(catname[0])
	all_catnames = [a_catnames, b_catnames, c_catnames]
	
else:
	all_catnames = [a_catnames, b_catnames]
	
	
	
	


for block_index,block in enumerate(all_catnames):


	#if a or b block, check there's 6 catalogue files (u, g, r, r2, i, NB)
	if block == all_catnames[0] or block == all_catnames[1]:
    		if len(block)!=6: 
    			import sys
    			print 'Block is incomplete: '
    			for line in block:
    				print line
    			#sys.exit()
    			#raw_input('Press any key to continue')



	for catpath in block:
	    	print "Processing ", catpath
    	
	   	#keep a copy of the original catalogue 
	  	orig_name = catpath[:-5]+'_original.fits'
	    	if not os.path.exists(orig_name):
	    		shutil.copyfile(catpath, orig_name)
	    		print 'Created ', orig_name   	
 	
    	
	    	#open the catalogue
	    	cat = fits.open(catpath, mode="update")
    	
    	   	
	    	#open corresponding ccd image to check the block/exposure filter
	        ccd_path, _ = catpath.rsplit('/', 1)
	        ccd_path, _ = ccd_path.rsplit('/', 1)
	        ccd_fpath = glob.glob(ccd_path+'/single/*decom.fit')
	        if len(ccd_fpath)==0:
	        	print 'Image fits file could not be found'
	        	print ccd_path
	        	import sys
	        	sys.exit()
	        
	        ccd = fits.open(ccd_fpath[0])   
          
	        #find out the filter used for the currect catalogue
	        flatsrc = ccd[1].header['FLATSRC']
	        filtername, _ = flatsrc.split('_', 1)
	        if filtername == 'r':
	            block_header = ccd[0].header['HIERARCH ESO OBS NAME']
	            _, block_header = block_header.rsplit('_', 1)
	            if block_header[0:2] == 'hr':
	                filtername = 'r'
	            elif block_header[0:2] == 'ur':
	                filtername = 'r2'
	        ccd.close()



    		#add useful columns to the catalogue: filtername, directory+ccd, exposure_time, magzpt, magnitude and error
    		for ccdnum in range(1,33):
        		print "ccd %i " %ccdnum
            
           		#delete the 'Blank' columns
   	   		#Could be better, but saves changing the rest of the script
   	   		#DO NOT USE: makes fits files invalid
           		#cols = cat[ccdnum].columns
           		#blank_names = [line for line in cols.names if 'Blank' in line]
   	   	 	#for name in blank_names:
   			#	cols.del_col(name)
   			#	cat[ccdnum].columns = cols

   		

        	    	table = cat[ccdnum].data
        	    	hdr = cat[ccdnum].header
        	    	colnames = table.dtype.names


        	   	if 'Filter' not in colnames:
        	        	#print "Adding filtername"
        	        	filtername_list = [filtername for x in range( len(table) )]
        	        	table = make_lists.append_table(table, 'Filter', filtername_list, np.dtype('S10'))
 
   		    	#filepath
        	    	if 'Fname' not in colnames:
        	    	   #print "Adding directory name"
        	    	    img_name, _ = catpath.rsplit('/', 1)
        	    	    img_name, _ = img_name.rsplit('/', 1)
        	    	    _ ,img_name = img_name.rsplit('/', 1)
        	    	    img_name_list = [img_name for x in range( len(table) )]
        	    	    table = make_lists.append_table(table, 'Fname', img_name_list, np.dtype('S20')) 


        	    	if 'Ccd_num' not in colnames:
        	    	    #print "Adding ccd number"
        	    	    num_list = [str(ccdnum) for x in range( len(table) )]
        	    	    table = make_lists.append_table(table, 'Ccd_num', num_list, np.dtype('>f4'))
        	           

			#APASSZPT = apass zpt (mags) for default extinction 
			#NIGHTZPT is the average apass zpt for the night
			if zptname not in colnames:   
      				zpt = hdr[zptname]         
      		        	magzpt_list = np.full( len(table) , zpt)
        	        	#print "Adding magnitude zeropoint"
        	        	table = make_lists.append_table(table, zptname, magzpt_list, np.dtype('>f4'))  

	
		  	#Adding exposure time column
		  	if 'Exp_time' not in colnames:
		    		expt = hdr['exptime']
        	        	expt_list = [expt for x in range( len(table) )]
        	        	#print "Adding exposure time"
        	        	table = make_lists.append_table(table, 'Exp_time', expt_list, np.dtype('>f4'))
         
                  

	    
        	    	#adding vega and AB magnitudes for each aperture
        	    	aperture_list = [name for name in table.dtype.names if 'Aper' in name and '_err' not in name and 'mag' not in name and 'corr' not in name and 'lower' not in name and 'upper' not in name]
        	    	for apname in aperture_list:
        	    
        	    		#skip apertures larger than 7
        	    		_, apnum = apname.rsplit('flux_')
        	    		try:
        	    			int(apnum)
        	     			if int(apnum) > 7 or int(apnum)<2:
        	    				continue
        	    		except: #_upper_lim
        	    			apnum, _ = apnum.split('_', 1)
        	    			apnum = int(apnum)
        	    			print 'Aperture', apnum
        	    			

            	          	
        	    		#uncorrected AB magnitudes
        	   	    	#print 'Adding ', apname, 'AB magnitide'	
        	    		AB_mags = [( -2.5*math.log10(line / hdr['exptime']) ) + hdr[zptname] if line>0 else float('nan')  for line in table[apname]  ]
        	    		if apname+'_mag_AB' not in colnames:
        	  		 	table = make_lists.append_table(table, apname+'_mag_AB', AB_mags, np.dtype('>f4'))
        	  		else:
        	  		 	table[apname+'_mag_AB']=AB_mags


        	    		#Vega magnitudes	        	
        	    		#print 'Adding vega magnitude'
        	    		#convert the AB magnitudes using the calculated conversions
		    		if filtername == 'r' or filtername == 'r2': conv = r_conv
		    		if filtername == 'i': conv = i_conv
		    		if filtername == 'u': conv = u_conv
		    		if filtername == 'g': conv = g_conv
		    		if filtername == 'NB':
		    			if ccdnum in range(0,9): #A
		   		        	conv = NBa_conv
		    		   	elif ccdnum in range(9,17): #D
		   				conv = NBd_conv
		        		elif ccdnum in range(17, 25):#B
		        			conv = NBb_conv
		        		elif ccdnum in range(25, 33): #C
		        			conv = NBc_conv
	                        	
  
		   		vega_mags = [line-conv for line in AB_mags]
		   		if apname+'_mag' not in colnames:
		        		table = make_lists.append_table(table, apname+'_mag', vega_mags, np.dtype('>f4'))  	
		        	else:
		        		table[apname+'_mag']=vega_mags


				#aperture corrected AB mags for every filter except u and NB
				if filtername != 'u' and filtername!='NB':
	
					#find the aperture correction txt file created by calc_ap_offset.py
        				if block_index==0: 
		        	   		corrections_fpath = os.getcwd()+'/aperture_corrections/a_'+filtername+'_aper'+str(apnum)+'_corrections.txt'
	        	   	    	elif block_index==1: 
	        	       			corrections_fpath = os.getcwd()+'/aperture_corrections/b_'+filtername+'_aper'+str(apnum)+'_corrections.txt'	        
	        	       		elif block_index==2: 
	        	       			corrections_fpath = os.getcwd()+'/aperture_corrections/c_'+filtername+'_aper'+str(apnum)+'_corrections.txt'       		
	
	        	     		if not os.path.exists(corrections_fpath):
	        	     			print 'ERROR: No aperture correction file'
	        	     			print corrections_fpath
	        	     			#raw_input('Press any key to continue')
	        	     			continue
	             			
	             			
        	       			#print 'Adding aperture corrected', apname, 'AB magnitide'
        	       			#read in the aperture corrections
               				corrections_list = []
            				with open(corrections_fpath) as f:
              					for line in f:
                					line=line.split()
                					corrections_list.append(float(line[1]))
                				
                			#calculate aperture corrected magnitudes			
                			if np.isnan(corrections_list[ccdnum-1]): 
                				corr_AB_mags = [ float('nan') for line in AB_mags]
                				corr_vega_mags = [ float('nan') for line in AB_mags]
                			else:		
                				#true_mag = line - apcor		
                				corr_AB_mags = [line-corrections_list[ccdnum-1] for line in AB_mags]  #0 counting	
                				corr_vega_mags = [line-conv for line in corr_AB_mags]  #0 counting
                			
                			
                			
                			
                			#append the table
                			if apname+'_corr_AB' not in colnames:
            					table = make_lists.append_table(table, apname+'_corr_AB', corr_AB_mags, np.dtype('>f4'))
            				else:
            					table[ apname+'_corr_AB' ] = corr_AB_mags
            					
 
                			if apname+'_corr' not in colnames:		
                				table = make_lists.append_table(table, apname+'_corr', corr_vega_mags, np.dtype('>f4'))
                 			else:
            					table[ apname+'_corr' ] = corr_vega_mags
                			        

                			     		
                			     		

                		     	
                	   

                                       	
           
			#save changes to the table 	
			cat[ccdnum].data = table
        	cat.flush()

    	cat.close(output_verify='ignore')
   

end = datetime.now()
print "Time taken: ", end-start

            
