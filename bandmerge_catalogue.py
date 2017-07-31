#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Bandmerges one vphas pointing into three catalogues using the a,b,c offsets. Removes blank cols"""

from astropy.io import fits
import glob
import os
import shutil
import operator
from datetime import datetime
import make_lists



#path to jystilts
jy_path = '/mirror/scratch/hbarker/pkgs/jystilts.jar'
#path to tmatchn.py. Is called as part of os, so must be the whole path
match_fpath = '/home/hbarker/scripts/tmatchn.py'
match2_fpath = '/home/hbarker/scripts/tmatch2.py'



args = make_lists.get_vphas_num()
ex_path = os.getcwd() + '/vphas_' + args.vphas_num + '_ex'
all_subdirs = [dirname for dirname in glob.glob(ex_path+'/*') if os.path.isdir(dirname)]
match_radius = 0.45 #arcsec

		
#Need to split this list into a,b,c then merge each set of catalogues.
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)

a_catnames = []
for dirname in a_block:
	catname = glob.glob(dirname + '/catalogues/*cat.fits')
	if len(catname)!=0:
		a_catnames.append(catname[0])


b_catnames = []
for dirname in b_block:
	catname = glob.glob(dirname + '/catalogues/*cat.fits')
	if len(catname)!=0:
		b_catnames.append(catname[0])

c_catnames = []
for dirname in c_block:

	if str(dirname) == 'None': 
		catname = 'None'
	else:
		catname = glob.glob(dirname + '/catalogues/*cat.fits')
		catname = catname[0]

	c_catnames.append(catname)




#match using tmatchn.py. Cols are in order u,g,r,r2,i,NB
#Had to redo into merging ccds with pair matching as couldn't fix max recusion limit with sys.setrecursionlimit() 
#Catalogue with the most entries is used as the reference in the pair matching
#Group matching now works: forgot to convert radians to degrees
savename = 'a_block_merged_cat.fits'
if not os.path.exists(os.getcwd()+'/'+savename):
	start = datetime.now()
	print 'Merging block a '

	#resave ccds as seperate files, pipe into tmatchn, delete files
	a_merged_dirname = os.getcwd() + '/a_merged_cat'
	if not os.path.exists(a_merged_dirname):
		os.makedirs(a_merged_dirname)
	
	for i in range(1,33): #loop over ccds
		print 'ccd ', str(i)
		
		savename = a_merged_dirname+'/ccd'+str(i)+'.fits'
		#if os.path.exists(savename):
			#continue
			#os.remove(savename)
		
		
		tempfpaths = []
		cat_lengths = []

		for j in range(0,6): #loop over the 6 filters and write ccds to file
			openfile = fits.open(a_catnames[j])
			hdr = openfile[i].header
			tab = openfile[i].data
			
			tempfpath = a_catnames[j][:-5]+'_ccd'+str(i)+'.fits'
			if os.path.exists(tempfpath):
				os.remove(tempfpath)
				
			tempfpaths.append(tempfpath)
			cat_lengths.append(len(tab))
			ccdfile = fits.BinTableHDU(header=hdr, data=tab)
			ccdfile.writeto(tempfpath)
		index, _ = max(enumerate(cat_lengths), key=operator.itemgetter(1)) #choose longest catalogue
		index+=1
		#merge catalogues and remove blankcols
		
		command_line = 'java -jar '+jy_path+' '+match_fpath+' {0} {1} {2} {3} {4} {5} {6} {7} {8}'
		os.system(command_line.format(tempfpaths[0], tempfpaths[1], tempfpaths[2], tempfpaths[3], tempfpaths[4], tempfpaths[5], index, match_radius, savename))

		
		#delete single filter ccd catalogues
		for fpath in tempfpaths:
			os.remove(fpath)
	end=datetime.now()
	print end-start

		
#merge the seperate ccd merged catalogues
if not os.path.exists(os.getcwd()+'/'+savename):
	print 'Merging a ccd catalogues'
	
	catnames = glob.glob(os.getcwd() + '/a_merged_cat/*.fits')
	merged_fname = os.getcwd()+'/a_block_merged_cat.fits'
	
	for ind, fname in enumerate(catnames):
		print 'ccd ', str(ind+1)
		
		if ind==0: 
			shutil.copyfile(fname, merged_fname)
		else:
			t1 = fits.open(merged_fname)
			t2 = fits.open(fname)
			
			nrows1 = t1[1].data.shape[0]
			nrows2 = t2[1].data.shape[0]
			
			nrows = nrows1+nrows2
			
			hdu = fits.BinTableHDU.from_columns(t1[1].columns, nrows=nrows)
			for colname in t1[1].columns.names:
				hdu.data[colname][nrows1:] = t2[1].data[colname]
			hdu.writeto(merged_fname, clobber=True)
			
			t1.close()
			t2.close()


print
savename = 'b_block_merged_cat.fits'
if not os.path.exists(os.getcwd()+'/'+savename):
	print 'Merging block b'	
	
	#resave ccds as seperate files, pipe into tmatchn, delete files
	b_merged_dirname = os.getcwd() + '/b_merged_cat'
	if not os.path.exists(b_merged_dirname):
		os.makedirs(b_merged_dirname)
	
	for i in range(1,33): #loop over ccds
		print 'ccd ', str(i)
		tempfpaths = []
		cat_lengths = []
		
		for j in range(0,6): #loop over the 6 filters
			openfile = fits.open(b_catnames[j])
			hdr = openfile[i].header
			tab = openfile[i].data
			
			tempfpath = b_catnames[j][:-5]+'_ccd'+str(i)+'.fits'
			if os.path.exists(tempfpath):
				os.remove(tempfpath)
				
			tempfpaths.append(tempfpath)
			cat_lengths.append(len(tab))
			
			ccdfile = fits.BinTableHDU(header=hdr, data=tab)
			ccdfile.writeto(tempfpath)
			
		savename = b_merged_dirname+'/ccd'+str(i)+'.fits'
		if os.path.exists(savename):
			os.remove(savename)
			
		index, _ = max(enumerate(cat_lengths), key=operator.itemgetter(1)) #choose longest catalogue
		index+=1
		
		command_line = 'java -jar '+jy_path+' '+match_fpath+' {0} {1} {2} {3} {4} {5} {6} {7} {8}'
		os.system(command_line.format(tempfpaths[0], tempfpaths[1], tempfpaths[2], tempfpaths[3], tempfpaths[4], tempfpaths[5], index, match_radius, savename))
		#remove single filter ccd catalogues 
		for fpath in tempfpaths:
			os.remove(fpath)

		print
	end=datetime.now()




if not os.path.exists(os.getcwd()+'/'+savename):
	print 'Merging b ccd catalogues'
	
	catnames = glob.glob(os.getcwd() + '/b_merged_cat/*.fits')
	merged_fname = os.getcwd()+'/b_block_merged_cat.fits'
	
	for ind, fname in enumerate(catnames):
		print 'ccd ', str(ind+1)
		
		if ind==0: 
			shutil.copyfile(fname, merged_fname)
		else:
			t1 = fits.open(merged_fname)
			t2 = fits.open(fname)
			
			nrows1 = t1[1].data.shape[0]
			nrows2 = t2[1].data.shape[0]
			
			nrows = nrows1+nrows2
			hdu = fits.BinTableHDU.from_columns(t1[1].columns, nrows=nrows)
			
			for colname in t1[1].columns.names:
				hdu.data[colname][nrows1:] = t2[1].data[colname]
			hdu.writeto(merged_fname, clobber=True)
			
			t1.close()
			t2.close()
	
	
	

savename = 'c_block_merged_cat.fits'
if len(c_block)==2 and not os.path.exists(os.getcwd()+'/'+savename): 

	if 'None' not in c_catnames:

		print 'Merging block c'	
		table1fpath = str(c_catnames[0])
		table2fpath = str(c_catnames[1])
		
		command_line = 'java -jar '+jy_path+' '+match2_fpath+'  {} {} {} {} {} {} {} {} {} {} {} {} {} {}'
		
		os.system(command_line.format(table1fpath, table2fpath, 'fits', 'True', 'fits', 'True', 'RA', 'DEC', 'RA', 'DEC', 0.3, '1and2', 'best', savename) ) 
print	
	
	
	
