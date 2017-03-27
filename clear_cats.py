#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#delete any catalogue files ina  directory that aren't originals

import os
import glob
import shutil

vphas_num = raw_input('vphas number: ')
ex_path = os.getcwd()+'/vphas_'+vphas_num+'_ex'

allcats = glob.glob(ex_path+'/*/catalogues/*.fits')

cats = [line for line in allcats if 'orig' not in line]

for line in cats:
	os.remove(line)
	
allcats = glob.glob(ex_path+'/*/catalogues/*.fits')

for line in allcats:
	newname, _ = line.rsplit('_original', 1)
	newname+='.fits'
	shutil.copyfile(line, newname)
	os.remove(line)
	
print 'Done'
