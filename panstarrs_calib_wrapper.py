#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""A wrapper to call all the scripts needed to calibrate vphas data"""


import os
import glob
import subprocess as sp
import shutil



#request number is the number supplied by the CASU website
reqNo = raw_input('Enter request number: ')

#the pointing number is used to label the directories
vphas_num = raw_input('Enter the vphas pointing number: ')




print "Downloading..."
sp.call(["download_vphas.sh", reqNo, vphas_num])
print
print "Sorting files..."
os.system("sort_vphas.py -v %s" %vphas_num)
print
print 'Please check the sorted files'
raw_input('Press any key to continue')
print
print "Decompressing files..."
os.system("decompress_vphas.py -v %s" %vphas_num)
print
print "Extracting ccds..."
os.system("extract_ccds.py -v %s" %vphas_num)	
print 
print 'Downloading PanSTARRS data'
os.system("get_panstarrs.py")
print
print 'Calculating aperture corrections with Pan-STARRS'
os.system("calc_panstarrs_ap_offset.py -v %s -b a" %vphas_num )
os.system("calc_panstarrs_ap_offset.py -v %s -b b" %vphas_num )
os.system("calc_panstarrs_ap_offset.py -v %s -b c" %vphas_num )
print

print
print 'Processing catalogues'
os.system("process_catalogues.py -v %s" %vphas_num)
print
print 'Bandmerging'
os.system("bandmerge_catalogue.py -v %s" %vphas_num)
print
print 'Correcting u band magnitudes'
os.system("u_corrections_auto.py -v %s -b a" %vphas_num )	
os.system("u_corrections_auto.py -v %s -b b" %vphas_num )
print 
print 'Calculating errors'
os.system("calc_errs.py -v %s -b a" %vphas_num )	
os.system("calc_errs.py -v %s -b b" %vphas_num )
os.system("calc_errs.py -v %s -b c" %vphas_num )
print
print

#delete the old bandmerged catalogues and remake them using the catalogues with errors
bandmerge_names = glob.glob(os.getcwd()+'/*_block_merged_cat.fits')
for name in bandmerge_names:
	os.remove(name)
dirnames = glob.glob('*/')
for name in dirnames:	
	if 'merged_cat' in name:
		shutil.rmtree(name)

	
print 'Bandmerging'
os.system("bandmerge_catalogue.py -v %s" %vphas_num)



#remove uneccessary directories
#dirnames = glob.glob('*/')
#for name in dirnames:
#	
#	if '_fitsfiles' in name:
#		shutil.rmtree(name)
#		
#	if 'merged_cat' in name:
#		shutil.rmtree(name)
		
	

		



print '-------------------------------COMPLETE-------------------------------------'


























	
