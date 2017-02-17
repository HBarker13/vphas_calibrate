#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python


#Extract extensions from decompressed vphas image and confidence FITS files.  Pulls out 32 ccds and saves them in a new directory.
#comment out one option

import os
from astropy.io import fits
import glob
import numpy as np
import sys
import argparse

def make_dir(dirname):
   if not os.path.exists(dirname):
       os.makedirs(dirname)
       print "Directory created: %s" %(dirname)
       print

import make_lists

args = make_lists.get_vphas_num()
current_path = os.getcwd()
ex_path = current_path +  '/vphas_'+args.vphas_num+'_ex'

#find decompressed fits files and created dirs for ccds
decom_files = glob.glob(ex_path+'/*/*/*decom*.fit')
for fpath in decom_files:
   ccd_dirpath, _ = fpath.rsplit('_', 1)
   _, fname = ccd_dirpath.rsplit('/', 1)
   ccd_dirpath = ccd_dirpath + '_ccds'
   if not os.path.exists(ccd_dirpath):
      os.makedirs(ccd_dirpath)
      print "Directory created: %s " %ccd_dirpath

#extract 32 ccds from each decom fits file and create file if it
#doesn't already exist
   open_decom = fits.open(fpath)
   for x in range(1,33):
      new_ccd = open_decom[x]
      newName = fname + '_ccd' + str(x) + '.fit'
      newPath = ccd_dirpath + '/' + newName
      if not os.path.exists(newPath):
         new_ccd.writeto(newPath)
         print "%s extracted" %(newName)
      else:
         print "%s already extracted" %newName
   open_decom.close()

   
