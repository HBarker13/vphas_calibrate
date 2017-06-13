#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#Script to organise downloaded vphas files. If an error is thrown when reading in the filelist, check the filelist for repeated/errenous entries

import numpy as np
import os
import shutil
import sys
from astropy.io import fits

import make_lists

def make_dir_tree(dirname):

    if not os.path.exists(dirname):
        os.makedirs(dirname)
        print "Created", dirname
        
    if not os.path.exists(dirname + '/single'):
        os.makedirs(dirname + '/single')
        print "Created", dirname, '/single'
        
    if not os.path.exists(dirname + '/calib'):
        os.makedirs(dirname + '/calib')
        print "Created", dirname,'/calib'
        
    if not os.path.exists(dirname + '/catalogues'):
        os.makedirs(dirname + '/catalogues')
        print "Created", dirname,'/catalogues' 






args = make_lists.get_vphas_num()
working_path  = os.getcwd()
vphas_dir =working_path + '/vphas_' + args.vphas_num 



#create a directory for the files to be sorted into
vphas_ex_path = working_path + '/vphas_' + args.vphas_num + '_ex'
if not os.path.exists( vphas_ex_path):
	os.makedirs( vphas_ex_path)
	print 'Created', vphas_ex_path




#read filelist from file into an array
filelist_path = vphas_dir + '/filelist_' + args.vphas_num

if not os.path.exists(filelist_path):
    print "Filelist cannot be found: ", filelist_path

with open(filelist_path) as fpath:
    filelist = [line.strip().split()[0] for line in fpath]
print "File list read in"



#Break up file list array into A, B and C, blocks.

#the filelist should have filenames in the order:
#single/xxxx
#catalogues/xxxx
#calib/xxxx
#calib/xxxx
#so get lines in groups of 4	


for x in range(0, len(filelist), 4):

	filegroup = filelist[x: x+4]

	#check the filenames to make sure they're in the expected order
	if 'single' not in filegroup[0]:
		print 'Filelist is not in the expected order'
		for l in filegroup: print l
		sys.exit()
	if 'catalogues' not in filegroup[1]:
		print 'Filelist is not in the expected order'
		for l in filegroup: print l
		sys.exit()
	if 'calib' not in filegroup[2]:
		print 'Filelist is not in the expected order'
		for l in filegroup: print l
		sys.exit()
	if 'calib' not in filegroup[3]:
		print 'Filelist is not in the expected order'
		for l in filegroup: print l
		sys.exit()


	#get the "root" filename. eg. o20240113_00120
	root = filegroup[0].lstrip('single/o').rstrip('.fit')

	#open the single file to see if it is in A, B, or C block
	#and get the filtername
	open_single = fits.open( vphas_dir + '/' + filegroup[0] )
	block_num = open_single[0].header['HIERARCH ESO TPL EXPNO']
	
	
	if block_num==1: block = 'a'
	elif block_num==2: block = 'b'
	elif block_num==3: block = 'c'
	
	
	#from Geert's email
	#hh => red concat, h-alpha
	#hr => red concat, r
	#hi => red concat, i

	#uu => blue concat, u
	#ug => blue concat, g
	#ur => blue concat, r
	
	
	filterlabel = open_single[0].header['HIERARCH ESO OBS NAME']
	_, filterlabel = filterlabel.rsplit('_', 1)
	
	if filterlabel[:2] == 'hr': filtername = 'r'
	elif filterlabel[:2] == 'ur': filtername = 'r2'	
	elif filterlabel[:2] == 'uu': filtername = 'u'
	elif filterlabel[:2] == 'ug': filtername = 'g'
	elif filterlabel[:2] == 'g': filtername = 'g'
	elif filterlabel[:2] == 'hi': filtername = 'i'
	elif filterlabel[:2] == 'hh': filtername = 'NB'
	else:
		print 'Error finding filter name'
		print filegroup[0]
		print open_single[0].header['HIERARCH ESO OBS NAME']
		raw_input('Press any key to continue')
		

	open_single.close()
	
	
	

	
	#create directories to sort this block into
	root_dir = vphas_ex_path + '/' + filtername + '_' + root + '_' + block
	make_dir_tree( root_dir)
	
	
	#move files into these directories
	for fname in filegroup:
	
		oldname = vphas_dir + '/' + fname
		newname = root_dir + '/' + fname
		
		if not os.path.exists(newname):
			shutil.copyfile(oldname, newname)
			print 'Copied', newname
	print

print
print "File sorting complete"

