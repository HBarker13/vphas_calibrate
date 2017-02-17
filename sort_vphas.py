#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#Script to organise downloaded vphas files. If an error is thrown when reading in the filelist, check the filelist for repeated/errenous entries

import numpy as np
import os
import shutil
from fnmatch import fnmatch
import sys
from astropy.io import fits

import make_lists

def make_dir_tree(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)
        print "Created %s" %(dirname)
    if not os.path.exists(dirname + '/single'):
        os.makedirs(dirname + '/single')
        print "Created %s" %(dirname + '/single')
    if not os.path.exists(dirname + '/calib'):
        os.makedirs(dirname + '/calib')
        print "Created %s" %(dirname + '/calib')
    if not os.path.exists(dirname + '/catalogues'):
        os.makedirs(dirname + '/catalogues')
        print "Created %s" %(dirname + '/catalogues') 

args = make_lists.get_vphas_num()

current_path  = os.getcwd()

#read filelist from file into an array
filelist_array = []
filelistPath = current_path + '/vphas_' + args.vphas_num + '/filelist_' + args.vphas_num

if not os.path.exists(filelistPath):
    print "Filelist cannot be found: %s" %filelistPath

with open(filelistPath) as fpath:
    for line in fpath:
        filelist_array.append(line[:-1]) #remove /n
print "File list read"


#Break up file list array into blocks. Assumes filelist is in order and there is a single, catalogue and two calib fields.
blockarray = []
for x in xrange(0, len(filelist_array), 4):
    block = filelist_array[x:x+4]
    code = block[0].lstrip('single/o').rstrip('.fit')
    block.insert(0,code)
    #checks first element is single and shows error if not
    if 'single' in block[1]: 
	blockarray.append(block)
    else:
        print "Problem"
        print block
        print
        print "Checking previous block: %s" % blockarray[-1][0]
        check_path = current_path + '/vphas_' + args.vphas_num + '/' + blockarray[-1][0]
        print check_path
        check_file = fits.open(check_path)
        print "File opened"
        header = check_file[0].header
        print "Correct filter: "
        esofilter = header['HIERARCH ESO INS FILT1 NAME']
        print esofilter
        print "Change %s then rerun" %filelistPath
        sys.exit()
 
 
for block in blockarray:
    #add colour to block[0] using calib file in block[3]
    colour = block[3][6:][0]
    if colour=='N':colour='NB'

    newDirPath = current_path + '/vphas_' + args.vphas_num  + '_ex/'+colour+'_'+ block[0]
    print "Making new directory: %s " %newDirPath
    make_dir_tree(newDirPath)


    for i in range(1,len(block)):
        old_path = current_path + '/vphas_' + args.vphas_num  + '/' + block[i]
        print old_path
        new_path = current_path + '/vphas_' + args.vphas_num + '_ex/'+colour+'_' + block[0] + '/' + block[i]
        print new_path
        if not os.path.exists(new_path):
            shutil.copyfile(old_path, new_path)
            print "Copied"
 
 
print "File sorting complete"

