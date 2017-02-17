#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python


#python script to decompress sorted files using imcopy command

import subprocess
import os
import sys
import glob

import make_lists

args = make_lists.get_vphas_num()
current_path = os.getcwd()
ex_path = current_path + '/vphas_' + args.vphas_num + '_ex'

dirnames = [fpath for fpath in glob.glob(ex_path+'/*') if 'trimmed' not in fpath and 'Halpha' not in fpath and '.txt' not in fpath]
for dirname in dirnames:
	print dirname
	#decompress single files
	singlepath = dirname+'/single'
	singlefile = glob.glob(singlepath+'/*.fit')
	singlefile= [fname for fname in singlefile if 'decom' not in fname]
	singlefile=singlefile[0]
	print singlefile
	decom_path = singlefile[:-4] + '_decom.fit'
	if not os.path.isfile(decom_path):
                subprocess.call(["imcopy", singlefile, decom_path ])
                print "Decompressed" 
        else:
                print "Already exists"
                continue
	
	#decompress calib files
	calibpath = dirname+'/calib'
	confpath = glob.glob(calibpath+'/*conf*.fit')
	conffile = [fname for fname in confpath if 'decom' not in fname]
	conffile = conffile[0]
	decom_path = conffile[:-4] + '_decom.fit'
	if not os.path.isfile(decom_path):
                subprocess.call(["imcopy", conffile, decom_path ])
                print "Decompressed"
        else:
                print "Already exists"
                continue
	
	
	
	
	
	
	
	

"""
#read filelist from the file into an array
filelist_list = []
filelistPath = current_path+'/vphas_'+args.vphas_num+ '/filelist_'+args.vphas_num
with open(filelistPath) as list:
    for line in list:
        filelist_list.append(line[:-1]) #remove /n when reading in

#Assumes filelist is in order and there is a single, catalogue and two calib fiels
blocklist = []
blockarray = []
counter = 0
for line in filelist_list:
    if 'single' in line:
        code = line.lstrip('single/o').rstrip('.fit')
        blocklist.append(code)
        blocklist.append(line)
        counter = 1
    if 'calib' in line:
        blocklist.append(line)
        counter += 1
    if 'catalogue' in line:
        blocklist.append(line)
        counter += 1
    if counter == 4:
        blockarray.append(blocklist)
        blocklist = []

for block in blockarray:
    for i in range(1,len(block)):
        print
        print "Processing %s" %(block[i])
        filepath = current_path +'/vphas_' + args.vphas_num + '_ex/' + block[0] + '/' + block[i]
        #skip catalogue files as they don't need decompressing
        if 'catalogue' in filepath:
            print "'catalogue' in filename %s.  File skipped" %(filepath)
        else:
            decom_path = filepath[:-4] + '_decom.fit'
            if not os.path.isfile(decom_path):
                subprocess.call(["imcopy", filepath, decom_path ])
                print "%s decompressed" %(block[i])
            else:
                print "%s already exists" %(decom_path)
"""
 
