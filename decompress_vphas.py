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
	
	if len(singlefile)==0: continue
	
	singlefile = singlefile[0]
	print singlefile
	decom_path = singlefile[:-4] + '_decom.fit'
	if not os.path.isfile(decom_path):
                subprocess.call(["imcopy", singlefile, decom_path ])
                print "Decompressed" 
        # else:
        #        print "Already exists"

                
	
	#decompress calib files
	calibpath = dirname+'/calib'
	confpath = glob.glob(calibpath+'/*conf*.fit')
	conffile = [fname for fname in confpath if 'decom' not in fname]
	

	if len(conffile)==0: continue
	
	conffile = conffile[0]
	print conffile
	decom_path = conffile[:-4] + '_decom.fit'
	if not os.path.isfile(decom_path):
                subprocess.call(["imcopy", conffile, decom_path ])
                print "Decompressed"
        else:
                print "Already exists"
                continue
	
	
	
	

