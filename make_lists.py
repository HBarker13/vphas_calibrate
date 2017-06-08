#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python
#functions to make lists of files/dirs in vphas directories


import glob
import os
import argparse
from astropy.io import fits
from astropy import wcs
import sys
import numpy as np
import math



#------------------------------------------------------------------------------------------------------#

#get args.vphas_num from command line
#use as: args = make_lists.get_vphas_num()
def get_vphas_num():
    parser = argparse.ArgumentParser(description="Collect vphas mosaicing and processing inputs")
    parser.add_argument('-v','--vphas_num', help="vphas pointing number", required=True)
    return parser.parse_args()
    
#------------------------------------------------------------------------------------------------------#    

#get all args from command line
def get_args():
    parser = argparse.ArgumentParser(description="Collect vphas mosaicing and processing inputs")
    parser.add_argument('-v','--vphas_num', help="vphas pointing number", required=True)
    parser.add_argument('-b','--bin_choice', help="bin pixels? y/n", required=True)
    parser.add_argument('-s','--bin_level', help="bin shrink factor (optional)", required=False)
    return parser.parse_args()

#------------------------------------------------------------------------------------------------------#

#list of paths to single, calib, cat dirs. Skip Halpha_div_R. Doesn't list fit files
#example output: [[expath/blah.png, expath/
def list_dirs(expath):
    dirlist = glob.glob(expath + '/*')
    paths_tree = []
    for dirpath in dirlist:           
        paths_branch = []
        paths_branch.append(dirpath)
        for subdir in glob.glob(dirpath + '/*'):
            paths_branch.append(subdir)
            for ccd_dir in glob.glob(subdir + '/*'):
                if not '.fit' in ccd_dir:
                    paths_branch.append(ccd_dir)
        paths_tree.append(paths_branch)
    return paths_tree

#------------------------------------------------------------------------------------------------------#

#lists colours available for mosaicking in expath (not Halpha/R, Halpaha-R)
def list_colours(expath):
    available_colours = set()
    dirnames = [x[1] for x in os.walk(expath)]
    topnames=dirnames[0]
    letters = [line[0:2] for line in topnames if line[0:2]!='tr' and line[0:2]!='Ha']
    for letter in letters:
    	available_colours.add(letter.strip('_'))
    return available_colours

#------------------------------------------------------------------------------------------------------#

#list all colours available for mosaicking in expath inc (Halpha_div_R...)
def list_colours_full(expath, bin_choice, bin_level):
    available_colours = set()
    dirnames = [x[1] for x in os.walk(expath)]
    topnames=dirnames[0]
    letters = [line[0:2] for line in topnames if line[0:2]!='tr' and line[0:2]!='Ha']
    for letter in letters:
    	available_colours.add(letter.strip('_'))

    Halpha_div_path = expath + '/Halpha_div_R'
    Halpha_sub_path = expath + '/Halpha_sub_R'
    if bin_choice == 'y':
        Halpha_div_path = Halpha_div_path + '/bin'+str(bin_level)
        Halpha_sub_path = Halpha_sub_path + '/bin'+str(bin_level)
    elif bin_choice == 'n':
        Halpha_div_path = Halpha_div_path + '/no_bin'
        Halpha_sub_path = Halpha_sub_path + '/no_bin'
    if os.path.exists(Halpha_div_path):
        available_colours.add('Halpha_div_R')
    if os.path.exists(Halpha_sub_path):
        available_colours.add('Halpha_sub_R')            
    return available_colours
 
#------------------------------------------------------------------------------------------------------#   

#make lists to branches of chosen colour (ugriNB). paths_tree from list_dirs
#eg. returns [root, root/confcorr, root/single, root/bin ....]
def colour_branches_wanted(paths_tree, colour_choice):
    if '_' not in colour_choice: colour_choice +='_'
    branches_to_process = []
    for branch in paths_tree:
        for leaf in branch:
            if 'calib' in leaf and colour_choice in leaf:
                if branch not in branches_to_process:
                    branches_to_process.append(branch)
    return branches_to_process

#------------------------------------------------------------------------------------------------------#
    
#list of top directories containing single, conf, confcorr etc. of a desired colour. NB. Colour must enter as 'r_', not 'r'
#eg. returns ~/vphas/vphas_01_ex/20140101_0001
def chosen_colour_pathlist(expath, colourChoice):
    if '_' not in colourChoice: colourChoice+='_'
    colour_list = set()
    if colourChoice == 'Halpha_div_R_':
        colour_list.add(expath+'/Halpha_div_R')
    if colourChoice == 'Halpha_sub_R_':
        colour_list.add(expath+'/Halpha_sub_R')
    else:
   	dirlist = glob.glob(expath + '/*') 
   	for dirpath in dirlist:          
            calib_path = dirpath + '/calib'
            calib_files = glob.glob(calib_path + '/*.fit')
            for fname in calib_files:
                if colourChoice in fname:
                    colour_list.add(dirpath)
    return colour_list

#------------------------------------------------------------------------------------------------------#

#input list of red and nb dirs that need pairing. Coord details for each are printed 
def print_dir_coords(colour1, colour1_dirlist):
        print "%s ccd1 CRPIX1.comment, CRPIX2.comment, CRPIX1, CRPIX2 values:" %(colour1)
        for dirname in colour1_dirlist:
                CCD = glob.glob(dirname+'/confcorr/*ccd1_*')
                openCCD = fits.open(CCD[0])
                refPix1 = openCCD[0].header['CRVAL1']
                refPix2 =  openCCD[0].header['CRVAL2']
                refPix1comment = [s for s in openCCD[0].header.comments['CRVAL1'].split()\
                                     if not s.isalpha() and s!='[deg]']
                refPix2comment = [s for s in openCCD[0].header.comments['CRVAL2'].split()\
                                     if not s.isalpha() and s!='[deg]']
                print dirname, refPix1comment, refPix2comment, refPix1, refPix2
                openCCD.close()
        print
        print


#------------------------------------------------------------------------------------------------------#

#add a column to a recarray
def append_table(table, name, arr, d_type):
    arr = np.asarray(arr)
    dtype = d_type
    newdtype = np.dtype(table.dtype.descr + [(name, d_type)])
    newtable = np.empty(table.shape, dtype=newdtype)
    for field in table.dtype.fields:
        newtable[field] = table[field]
    newtable[name] = arr
    return newtable


#-------------------------------------------------------------------------------------------------------#



#takes an uncollimated table and converts into recarray
#eg. tab = [[a[1], b[1], c[1]], [a[2], b[2], c[2]]    
#    r_array=[[a[1], a[2]], [b[1], b[2]], [c[1], c[2]] 
def make_recarray(tab, title_list):	
	dtype_list = ['>f4' for item in title_list]
	str_list = ['vphas_num', 'png', 'pn_name','name', 'v_name', 'PN', 'data_release', 'sourceID', 'primaryID', 'warning_u', 'detectionID_u', 'warning_g', 'detectionID_g', 'warning_r2', 'detectionID_r2', 'warning_ha', 'detectionID_ha', 'warning_r', 'detectionID_r', 'warning_i', 'detectionID_i', 'field', 'spectype', 'Companion_SpecType_(i)', 'Lower_SpecType_(i)', 'Upper_SpecType_(i)', 'Companion SpecType (J)', 'Lower SpecType (J)', 'Upper SpecType (J)', 'Abundance', 'filtername', 'Filter_1', 'Filter_2', 'Filter_3', 'Filter_4', 'Filter_5', 'Filter_6', 'Fname_1', 'Fname_2', 'Fname_3', 'Fname_4', 'Fname_5', 'Fname_6', 'pn', 'block', 'Lum_class']
	for ind, val in enumerate(title_list):
		if val in str_list:
			dtype_list[ind]='|S20'
			
	if str_override==True:
		dtype_list = ['|S20' for item in title_list]
			
	name_dtype = [tuple(line) for line in zip(title_list, dtype_list)]

	data_array = []
	for i in range(len(title_list)):
		col = [line[i] for line in tab]
		data_array.append(col)

	r_array = np.rec.fromarrays((data_array), dtype=name_dtype)
	return r_array	





#------------------------------------------------------------------------------------------------------#





#calculate vega magnitude from counts as mag = -2.5*log10(counts/exposure_time) + nightzpt
def calc_vega_mags(table, apname, appendix):
	mags = [(-2.5*math.log10(line[apname+'_'+appendix] / line['Exp_time_'+appendix]) ) + line['Magzpt_'+appendix] for line in table if line[apname+'_'+appendix]>0 ]
	return mags


#-------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------#

#calculate vega magnitude from counts as mag = -2.5*log10(counts/exposure_time) + nightzpt
def calc_vega_mags2(table, apname, exp_time, magzpt):
	mags = [(-2.5*math.log10(line[0] / line[1]) ) + line[2] for line in zip(table[apname], exp_time, magzpt) if line[0]>0 ]
	return mags


#-------------------------------------------------------------------------------------------------------#



#calculate AB mags from vega 
def calc_ab_mags(vega, filtername, ccdnum):
	u_conv = 2.5*math.log((3631.0/1550.81) ,10)
	g_conv = 2.5*math.log((3631.0/3960.53) ,10)
	r_conv = 2.5*math.log((3631.0/3094.68) ,10)
	i_conv = 2.5*math.log((3631.0/2563.84) ,10)
	NBa_conv = 2.5*math.log((3631.0/2903.53) ,10)
	NBb_conv = 2.5*math.log((3631.0/2931.14) ,10)
	NBc_conv = 2.5*math.log((3631.0/2929.13) ,10)
	NBd_conv = 2.5*math.log((3631.0/2669.32) ,10)
	
	if filtername == 'r_r': conv = r_conv
	if filtername == 'r_b': conv = r_conv
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
	ab_mags = [line+conv for line in vega]
	
	return ab_mags



#-------------------------------------------------------------------------------------------------------#

#save fits table
"""
opencat = fits.open(cat)
all_hdus = []
hdr = opencat[0].header
prihdu = fits.PrimaryHDU(header=hdr)
	all_hdus.append(prihdu)
	for i in range(1,33):
		print 'ccd ', str(i)
		
		#update table
		table = opencat[i].data
		colnames = table.dtype.names
		colfields = table.dtype.fields
		newtable = table[table['Classification']==-1]
		new_cols = []
		for name in colnames:
			if 'Blank' in name: continue
			fmt = colfields[name][0]
			arr = newtable[name]
			col = fits.Column(name=name, format=fmt, array=arr)
			new_cols.append(col)
				
		tbhdu = fits.BinTableHDU.from_columns(new_cols)
		all_hdus.append(tbhdu)
	tbhdulist = fits.HDUList(all_hdus)
	tbhdulist.writeto(newname, clobber=True)
	opencat.close()
"""
	
#-------------------------------------------------------------------------------------------------------#

#returns polynomial equation that can be used to calculate line of best fit
#uses only the section of the x axis between xmin and xmax
def calc_lobf(x, xmin, xmax, y):
	coords = [[line[0],line[1]] for line in zip(x,y) if xmin<x<xmax]
	newx = [line[0] for line in coords]
	newy = [line[1] for line in coords]
	coeffs = np.polyfit(newx, newy, 1)
	polynomial = np.poly1d(coeffs)
	return polynomial
	#bestfit = polynomial(x_axis_vals)
#-------------------------------------------------------------------------------------------------------#


#Pair directories overlapping ccds of different colours. Name directories with colours in order u,g,r,i,NB
#eg. matching u and r directories: colour1=u, colour2=r, or matching r and i directores: colour1=r, colour2=i
def make_pair_list(expath, colour1, colour2):
    dirlist = glob.glob(expath + '/*')
    keep_choice = 'n'
    #load pairlist.txt if it exists and read in 
    print "Searching for file containing %s and %s pair list..." %(colour1, colour2)
    pair_list = []
    pair_list_path = expath + '/'+colour1+'_'+colour2+'_pairlist.txt'
    if os.path.exists(pair_list_path):
        with open(pair_list_path) as pairfile:
            for line in pairfile:
                pair = []
                pair0, pair1 = line.split(',', 1)
                pair.append(pair0)
                pair.append(pair1[:-1])
                if pair not in pair_list:
                    pair_list.append(pair)
        pairfile.close()
        print "Pair list from file: "
        for line in pair_list:
                print line
        print 
        pair_choice = ''
        while pair_choice != 'y' and pair_choice != 'n':
                pair_choice = raw_input("Is this pair list correct? (y/n): ")
                print
        if pair_choice == 'y':
            return pair_list
            
    colour1Dirs = set()
    colour2Dirs = set()
    for dirpath in dirlist:
        calib_path = dirpath + '/calib'
        for fpath in glob.glob(calib_path+'/*'):
            _, fname = fpath.rsplit('/', 1)
            colour, _ = fname.split('_', 1)
            if colour == colour1:
                colour1Dirs.add(dirpath)
            if colour == colour2:
                colour2Dirs.add(dirpath)

    #file containing pair list either non-existent or rejected
    print "Searching for pairs using CRVAL header information..."
    pair_choice = ''
    pair_list = []
    print_dir_coords(colour1, colour1Dirs)
    print_dir_coords(colour2, colour2Dirs)
    if len(colour1Dirs) != len(colour2Dirs):
        print "WARNING: the number of %s dirs is not the same as %s dirs." %(colour1, colour2)
        print "%sDirs: %i" %(colour1, len(colour1Dirs))
        print "%sDirs: %i" %(colour2, len(colour2Dirs))

    #match directories of overlapping ccds using ra/dec in CRVAL comments
    for colour1Dir in colour1Dirs:
        colour1CCD = glob.glob(colour1Dir+'/confcorr/*ccd1_*') 
        openColour1 = fits.open(colour1CCD[0])
        colour1Commentpix1 = [s for s in openColour1[0].header.comments['CRVAL1'].split()\
                          if not s.isalpha() and s!='[deg]']
        colour1Commentpix2 = [s for s in openColour1[0].header.comments['CRVAL2'].split()\
                          if not s.isalpha() and s!='[deg]']
        colour1Valuepix1 = round(openColour1[0].header['CRVAL1'], 2)
        colour1Valuepix2 =  round(openColour1[0].header['CRVAL2'], 2)
        openColour1.close()
        for colour2Dir in colour2Dirs:
            colour2CCD = glob.glob(colour2Dir+'/confcorr/*ccd1_*') 
            openColour2 = fits.open(colour2CCD[0])
            colour2Commentpix1 = [s for s in openColour2[0].header.comments['CRVAL1'].split()\
                             if not s.isalpha() and s!='[deg]']
            colour2Commentpix2 = [s for s in openColour2[0].header.comments['CRVAL2'].split()\
                             if not s.isalpha() and s!='[deg]']
            colour2Valuepix1 = round(openColour2[0].header['CRVAL1'], 2)
            colour2Valuepix2 = round(openColour2[0].header['CRVAL2'], 2)
            openColour2.close()
            if colour2Commentpix1 == colour1Commentpix1 and colour2Valuepix1 == colour1Valuepix1:
                if colour2Commentpix2 == colour1Commentpix2 and colour2Valuepix2 == colour1Valuepix2 :
                    pair = [colour1Dir, colour2Dir]
                    #red directory has already been paired
                    if pair not in pair_list:
                        pair_list.append(pair)

    if len(pair_list) < len(colour1Dirs):
        print
        print "The number of pairs found is not equal to the number of %s directories." %(colour1)
        print
    print "Pair list generated by looking at CRVAL header entries: "
    for line in pair_list:
        print line
        print
    
    pair_choice = ''
    while pair_choice != 'y' and pair_choice != 'n':
        pair_choice = raw_input("Is this pair list correct? (y/n): ")
    if pair_choice == 'y':
        #add any dirs not in pair list to remainder_lists
        remainder_colour1dirs = set()
        remainder_colour2dirs = set()
        for dirname in colour1Dirs:
            check = 'n'
            for line in pair_list:
                if dirname == line[0]:
                    check = 'y'
            if check == 'n':
                remainder_colour1dirs.add(dirname)
        for dirname in colour2Dirs:
            check = 'n'
            for line in pair_list:
                if dirname == line[1]:
                    check = 'y'
            if check == 'n':
                remainder_colour2dirs.add(dirname)
        check_pair_list(expath, pair_list, remainder_colour1dirs, remainder_colour2dirs, pair_list_path)
        
        manual_choice = ''
        print
        while manual_choice != 'y' and manual_choice != 'n':
            manual_choice = raw_input("Continue and finish table manually (y) or save list of paired directories and exit (n)? ")
        if manual_choice == 'n':
                pairfile = open(pair_list_path, 'w')
                for line in pair_list:
                    line = str(line[0]) + ',' + str(line[1]) + '\n'
                    pairfile.write(line)
                pairfile.close
                print "Pair list saved "
                return pair_list


#elif pair_choice = 'n' or manual_choice = 'n'
    print "Table must be created manually"
    keep_saved_list = ''
    while keep_saved_list != 'y' and keep_saved_list != 'n':
        keep_saved_list = raw_input("Extend current list (y) or delete it (n)? ")
        if keep_saved_list == 'n':
            pair_list = []
    finish_choice = 'n'
    while finish_choice == 'n':
        print
        print_dir_coords(colour1, colour1Dirs)
        print_dir_coords(colour2, colour2Dirs)
        edit_choice = ''
        while edit_choice != 'q':
            edit_choice = ''
            numPairs = min(len(colour1Dirs), len(colour2Dirs))
            print
            print "Expecting %i pairs" %numPairs
            print "Current pair number: %i " %len(pair_list)
            print
            for line in pair_list:
                print line
            print
            while edit_choice != 'a' and edit_choice != 'd' and edit_choice != 'q':
                edit_choice = raw_input("Add pair (a), delete line (d) or end table creation (q)? ")

            if edit_choice == 'a':
                pair = []
                print "Enter pair: "
                pair0 = raw_input("Enter "+colour1+" Dir: " + expath+ '/')
                pair1 = raw_input("Enter "+colour2+" Dir: " + expath + '/' )
                pair.append(expath + '/' +pair0)
                pair.append(expath + '/' +pair1)
                pair_check = ''
                while pair_check != 'y' and pair_check != 'n':
                    pair_check = raw_input("Is this pair correct? (y/n) " + str(pair) + '  ')
                    if pair_check == 'y':
                        pair_list.append(pair)
                    if pair_check == 'n':
                        print "Pair will not be added to list"
                        print

            if edit_choice == 'd':
                colour1_to_del = expath +'/'+ raw_input("Enter the first directory in the pair to be deleted: "+expath+"/")
                colour2_to_del = expath + '/'+ raw_input("Enter the second directory in the pair to be deleted: "+expath+"/")
                try:
                    for pair in pair_list:
                        if pair[0] == colour1_to_del and pair[1] == colour2_to_del:
                            print pair
                            pair_list.remove(pair)
                            print "Pair removed"
                except:
                    print "Pair could not be found."


        #edit_choice =q so table creation ends and checking begins
        #add any dirs not in pair list to remainder_lists
        remainder_colour1dirs = set()
        remainder_colour2dirs = set()
        for dirname in colour1Dirs:
            check = 'n'
            for line in pair_list:
                if dirname == line[0]:
                    check = 'y'
            if check == 'n':
                remainder_colour1dirs.add(dirname)
        for dirname in colour2Dirs:
            check = ''
            for line in pair_list:
                if dirname == line[1]:
                    check = 'y'
            if check == 'n':
                remainder_colour2dirs.add(dirname)

        check_pair_list(expath, pair_list, remainder_colour1dirs, remainder_colour2dirs, pair_list_path)

        finish_choice = ''
        while finish_choice != 'y' and finish_choice != 'n':
            finish_choice = raw_input("Save pair list and exit table creation? (y/n) ")
        if finish_choice == 'y':
            print "Saving..."
            pairfile = open(pair_list_path, 'w')
            for line in pair_list:
                line = str(line[0]) + ',' + str(line[1]) + '\n'
                pairfile.write(line)
            pairfile.close
            return pair_list
            print
            print
        
#------------------------------------------------------------------------------------------------------------#
#check the pairs made in make pair list are correct. If colour2's ccd1 is more than 10 pixels in x
#or y from colour1's ccd1, the pairing is removed

def check_pair_list(expath, pair_list, remainder_colour1dirs, remainder_colour2dirs, savepath):
    #check shiftx, shifty for ccd1 in each directory
    print "Checking pair list to make sure shifts between ccds are <10 pixels..."
    for pair in pair_list:
        print
        print pair[0]
        print pair[1]
        colour1_ccd = glob.glob(pair[0]+'/confcorr/*_ccd1_*.fit')
        colour2_ccd = glob.glob(pair[1]+'/confcorr/*_ccd1_*.fit')

        colour1 = fits.open(colour1_ccd[0])
        colour2 = fits.open(colour2_ccd[0])
        colour1_header = colour1[0].header
        colour2_header = colour2[0].header
        colour1_wcs = wcs.WCS(colour1_header)
        colour2_wcs = wcs.WCS(colour2_header)
        colour1_img = colour1[0].data
        colour2_img = colour2[0].data

        cx = colour1_img.shape[0]/2
        cy = colour1_img.shape[1]/2
        cra, cdec = colour1_wcs.wcs_pix2world(cx, cy, 1) #colour1 as reference
        refx, refy = colour2_wcs.wcs_world2pix(cra, cdec, 1) #pixels in colour2 with same wcs as cra and cdec in colour1
        shiftx = refx-cx
        shifty = refy-cy
        print "shiftx, shifty: ", shiftx, shifty
        colour1.close()
        colour2.close()
            
        if abs(shiftx)>10 or abs(shifty)>10: 
            #shifts are too large pairing must be wrong
            print "Shift > 20 between: "
            print pair[0]
            print pair[1]
            print "Directories must be incorrectly paired. Will attempt to divide remaining pairs and return to these at the end."
            remainder_colour1dirs.add(pair[0])
            remainder_colour2dirs.add(pair[1])
            pair_list.remove(pair)

    if len(remainder_colour1dirs) == 0 and len(remainder_colour2dirs) == 0:
        print "No remainder directories. "
        print "Pair list has been created and saved: "
        pairfile = open(savepath, 'w')
        for line in pair_list:
            line = str(line[0]) + ',' + str(line[1]) + '\n'
            pairfile.write(line)
        pairfile.close
        for line in pair_list:
            print line
        return pair_list
        
    else:
        print
        print "Checking completed. Remainder unpaired directories exist. "
        return remainder_colour1dirs
        return remainder_colour2dirs
        print 
        print_dir_coords(remainder_colour1dirs, remainder_colour2dirs)
        print 
        

#-----------------------------------------------------------------------------------------------------#
#match directories in all 5 colours that contain ccds that overlap. Assumes u,g,r,i,Halpha all exist.
#returns list with colours in order r,nb,i,u,g

def match_all_colours(expath):
    dirlist = glob.glob(expath + '/*')
    keep_choice = 'n'
    #load pairlist.txt if it exists and read in 
    print "Searching for file containing matches of all colours..."
    match_list = []
    fpath = expath + '/all_colour_matchlist.txt'
    if os.path.exists(fpath):
        with open(fpath) as f:
            for line in f:
                matches = line.split()
                if matches not in match_list:
                    match_list.append(matches)
        f.close()
        print "List from file: "
        for line in match_list:
                print line
        print 
        pair_choice = ''
        while pair_choice != 'y' and pair_choice != 'n':
                pair_choice = raw_input("Is this match list correct? (y/n): ")
                print
        if pair_choice == 'y':
            return match_list
    
    u_dirs = chosen_colour_pathlist(expath, 'u_')
    r_dirs = chosen_colour_pathlist(expath, 'r_')
    nb_dirs = chosen_colour_pathlist(expath, 'NB_')
    i_dirs = chosen_colour_pathlist(expath, 'i_')
    g_dirs = chosen_colour_pathlist(expath, 'g_')
    colour_dirs_list = [nb_dirs, i_dirs, u_dirs, g_dirs]

    #file containing pair list either non-existent or rejected
    print "Searching for pairs using CRVAL header information..."
    match_choice = ''
    print
    print_dir_coords('r', r_dirs)
    print_dir_coords('NB', nb_dirs)
    print_dir_coords('i', i_dirs)
    print_dir_coords('g', g_dirs)
    print_dir_coords('u', u_dirs)

    #match directories of overlapping ccds using ra/dec in CRVAL comments using red as reference
    #less strict matching criteria than pair matching above
    match_list = []
    for r_dirname in r_dirs:
        match = []
        redCCD = glob.glob(r_dirname+'/confcorr/*ccd1_*') 
        red = fits.open(redCCD[0])
        rCommentpix1 = [s for s in red[0].header.comments['CRVAL1'].split()\
                          if not s.isalpha() and s!='[deg]']
        rCommentpix2 = [s for s in red[0].header.comments['CRVAL2'].split()\
                          if not s.isalpha() and s!='[deg]']
        rValuepix1 = round(red[0].header['CRVAL1'], 1)
        rValuepix2 =  round(red[0].header['CRVAL2'], 1)
        red.close()
        match.append(r_dirname)
        for colour in colour_dirs_list:
            for colour2Dir in colour:
                colour2CCD = glob.glob(colour2Dir+'/confcorr/*ccd1_*') 
                openColour2 = fits.open(colour2CCD[0])
                colour2Commentpix1 = [s for s in openColour2[0].header.comments['CRVAL1'].split()\
                                      if not s.isalpha() and s!='[deg]']
                colour2Commentpix2 = [s for s in openColour2[0].header.comments['CRVAL2'].split()\
                if not s.isalpha() and s!='[deg]']
                colour2Valuepix1 = round(openColour2[0].header['CRVAL1'], 1)
                colour2Valuepix2 = round(openColour2[0].header['CRVAL2'], 1)
                openColour2.close()
                if colour2Commentpix1 == rCommentpix1 or colour2Valuepix1 == rValuepix1:
                    if colour2Commentpix2 == rCommentpix2 or colour2Valuepix2 == rValuepix2 :
                        if colour2Dir not in match:
                            match.append(colour2Dir)
                            continue #makes sure only one directory of each colour is added
        if len(match) == 5: #all colours have been matched
            match_list.append(match)

    if len(match_list) == 0:
        print "All colours could not be matched"

    print "Pair list(s) generated by looking at CRVAL header entries: "
    for i in range(0, len(match_list)):
        num = str(i)+'. '
        print 
        print num, match_list[i]
    
    print
    pair_choice = ''
    while pair_choice != 'y' and pair_choice != 'n':
        pair_choice = raw_input("Check this list? (y/n): ")
    if pair_choice == 'y':
        checked_matches = []
        for i in range(len(match_list)):
            #add any dirs not in pair list to remainder_lists
            remainder_rdirs = set([dirname for dirname in r_dirs for line in match_list if dirname!=line[0]])
            remainder_nbdirs = set([dirname for dirname in nb_dirs for line in match_list if dirname!=line[1]])
            remainder_idirs = set([dirname for dirname in i_dirs for line in match_list if dirname!=line[2]])
            remainder_udirs = set([dirname for dirname in u_dirs for line in match_list if dirname!=line[3]])
            remainder_gdirs = set([dirname for dirname in g_dirs for line in match_list if dirname!=line[4]])
     
            print "Checking match ", i
            match = check_fivematch_list(expath, match_list[i], remainder_rdirs, remainder_nbdirs, remainder_idirs,\
                                              remainder_udirs, remainder_gdirs)
            if match != []:
                check_matches.append(match)
        match_list = checked_matches
                

        if len(match_list)==0:
            print "No matches are correct. "
        else:
            manual_choice = ''
            while manual_choice != 'y' and manual_choice != 'n':
                manual_choice = raw_input("Continue and finish table manually (y) or save list of paired directories and exit (n)? ")
            if manual_choice == 'n':
                savefile = open(fpath,  'w')
                for line in match_list:
                    new_line = ''
                    for i in range(0, len(line)):
                        new_line += str(line[i])
                        new_line += ' '
                    savefile.write(new_line)
                savefile.close
                print "Match list saved "
                return match_list


#elif pair_choice = 'n' or manual_choice = 'n'
    print "Table must be created manually"
    if len(match_list)>0:
        keep_saved_list = ''
        while keep_saved_list != 'y' and keep_saved_list != 'n':
            keep_saved_list = raw_input("Extend current list (y) or delete it (n)? ")
            if keep_saved_list == 'n':
                pair_list = []
    finish_choice = 'n'
    while finish_choice == 'n':
        print
        print_dir_coords('r', r_dirs)
        print_dir_coords('NB', nb_dirs)
        print_dir_coords('i', i_dirs)
        print_dir_coords('g', g_dirs)
        print_dir_coords('u', u_dirs)
        edit_choice = ''
        while edit_choice != 'q':
            edit_choice = ''
            for i in range(0, len(match_list)):
                print i, '. ', match_list[i]
            print
            while edit_choice != 'a' and edit_choice != 'd' and edit_choice != 'q':
                edit_choice = raw_input("Add match (a), delete line (d) or end table creation (q)? ")

            if edit_choice == 'a':
                match = []
                print "Enter directories: "
                r_match = raw_input("Enter r dir: " + expath+ '/')
                nb_match = raw_input("Enter nb dir: " + expath + '/' )
                i_match = raw_input("Enter i dir: " + expath+ '/')
                u_match = raw_input("Enter u dir: " + expath + '/' )
                g_match = raw_input("Enter g dir: " + expath+ '/')

                match.append(expath + '/' + r_match)
                match.append(expath + '/' + nb_match)
                match.append(expath + '/' + i_match)
                match.append(expath + '/' + u_match)
                match.append(expath + '/' + g_match)

                match_check = ''
                while match_check != 'y' and match_check != 'n':
                    match_check = raw_input("Is this correct? (y/n) " + str(match) + '  ')
                    if match_check == 'y':
                        match_list.append(match)
                    if match_check == 'n':
                        print "Will not be added to list"
                        print

            if edit_choice == 'd':
                match_to_del = int(raw_input("Enter number index of match to be deleted (starts at zero)."))
                del match_list[match_to_del]
                        

        #edit_choice =q so table creation ends and checking begins
        #add any dirs not in pair list to remainder_lists
        for match in match_list:
            remainder_rdirs = set([dirname for dirname in r_dirs for line in match_list if dirname!=line[0]])
            remainder_nbdirs = set([dirname for dirname in nb_dirs for line in match_list if dirname!=line[1]])
            remainder_idirs = set([dirname for dirname in i_dirs for line in match_list if dirname!=line[2]])
            remainder_udirs = set([dirname for dirname in u_dirs for line in match_list if dirname!=line[3]])
            remainder_gdirs = set([dirname for dirname in g_dirs for line in match_list if dirname!=line[4]])

            check_fivematch_list(expath, match,remainder_rdirs, remainder_nbdirs,\
                                 remainder_idirs, remainder_udirs, remainder_gdirs)

        finish_choice = ''
        while finish_choice != 'y' and finish_choice != 'n':
            finish_choice = raw_input("Save pair list and exit table creation? (y/n) ")
            if finish_choice == 'y':
                savefile = open(fpath,  'w')
                for line in match_list:
                    new_line = ''
                    for i in range(0, len(line)):
                        new_line += str(line[i])
                        new_line += ' '
                        savefile.write(new_line)
                savefile.close
                print "Match list saved "
                return match_list
                
#----------------------------------------------------------------------------------------------------#
#creates lists of the a,b and c pointings in one vphas directory. Returns in order u,g,r_r,r_b,i,NB
def bandmerge_list(ex_path):
	dirpaths = [dirname for dirname in glob.glob(ex_path+'/*') if os.path.isdir(dirname) and 'trimmed' not in dirname and 'Halpha' not in dirname]
	
	#split blocks into filternames, including red and blue r pointings
	i_dirs = []
	NB_dirs = []
	u_dirs = []
	g_dirs = []
	r_r_dirs = []
	r_b_dirs = []
	
	for dirpath in dirpaths:
		_, dirname = dirpath.rsplit('/', 1)
		filtername = dirname[0]
		if filtername == 'i' : i_dirs.append(dirpath)
		if filtername == 'u' : u_dirs.append(dirpath)
		if filtername == 'g' : g_dirs.append(dirpath)
		if filtername == 'N' : NB_dirs.append(dirpath)
		if filtername == 'r':
			blocktest_path = [fname for fname in glob.glob(dirpath+'/single/*.fit') if 'decom' not in fname]
			header = fits.getheader(blocktest_path[0])
			filterlabel = header['HIERARCH ESO OBS NAME']
			_, filterlabel = filterlabel.rsplit('_', 1)
			if filterlabel[:2] == 'hr': r_r_dirs.append(dirpath)
			if filterlabel[:2] == 'ur': r_b_dirs.append(dirpath)
			
			
	a_block = []
	b_block = []		
	#check the number of directories in each filter block
	#add to blocks, assume the lowest directory number in each filter is 'a', the highest is 'b' and the middle is 'c'
	#unless the name has been changed
	if len(u_dirs)!=2: 
		if len(u_dirs)==0:
			print 'No u directories: ', len(u_dirs)
		 	cont = raw_input('Press any key to continue')
		 	a_block.append('nan')
		 	b_block.append('nan')
		else:
			print 'Incorrect number of u directories: ', len(u_dirs)
			a_block.append(min(u_dirs))
			b_block.append(max(u_dirs))
	elif len(u_dirs)==2:
		a_block.append(min(u_dirs))
		b_block.append(max(u_dirs))
		
	if len(g_dirs)!=3: 
		if len(g_dirs)==0:
			print 'No g directories: ', len(g_dirs)
			cont = raw_input('Press any key to continue')
			a_block.append('nan')
			b_block.append('nan')
		else:
			print 'Incorrect number of g directories: ', len(g_dirs)
			a_block.append(min(g_dirs))
			b_block.append(max(g_dirs))
	elif len(g_dirs)==3:
		a_block.append(min(g_dirs))
		check=False
		for gname in g_dirs:
			if gname[-1]=='b':
				b_block.append(gname)
				check = True
		if check==False:
			b_block.append(max(g_dirs))
		
	
			
	if len(r_r_dirs)!=2: 
		if len(r_r_dirs)==0:
			print 'No r_r directories: ', len(r_r_dirs)
		 	cont = raw_input('Press any key to continue')
		 	a_block.append('nan')
		 	b_block.append('nan')
		else:
			print 'Incorrect number of r_r directories: ', len(r_r_dirs)
			a_block.append(min(r_r_dirs))
			b_block.append(max(r_r_dirs))
	elif len(r_r_dirs)==2:
		a_block.append(min(r_r_dirs))
		b_block.append(max(r_r_dirs))


	if len(r_b_dirs)!=2: 
		if len(r_b_dirs)==0:
			print 'No r_b directories: ', len(r_b_dirs)
		 	cont = raw_input('Press any key to continue')
		 	a_block.append('nan')
		 	b_block.append('nan')
		else:
			print 'Incorrect number of r_b directories: ', len(r_b_dirs)
			a_block.append(min(r_b_dirs))
			b_block.append(max(r_b_dirs))
	elif len(r_b_dirs)==2:
		a_block.append(min(r_b_dirs))
		b_block.append(max(r_b_dirs))
		
			
	if len(i_dirs)!=2: 
		if len(i_dirs)==0: 
			print 'No i directories: ', len(i_dirs)
			cont = raw_input('Press any key to continue')
			a_block.append('nan')
			b_block.append('nan')
		else:
			print 'Incorrect number of i directories: ', len(i_dirs)
			a_block.append(min(i_dirs))
			b_block.append(max(i_dirs))
	elif len(i_dirs)==2:
		a_block.append(min(i_dirs))
		b_block.append(max(i_dirs))
	
	
	if len(NB_dirs)!=3:
		if len(NB_dirs)==0:
			print 'No NB directories: ', len(NB_dirs)
			cont = raw_input('Press any key to continue')
			a_block.append('nan')
			b_block.append('nan')
		else:
			print 'Incorrect number of NB directories: ', len(NB_dirs)
			a_block.append(min(NB_dirs))
			b_block.append(max(NB_dirs))
	elif len(NB_dirs)==3:
		a_block.append(min(NB_dirs))
		b_block.append(max(NB_dirs))

		
	#create c_block
	c_block = []
	for dirname in NB_dirs:
		if dirname not in a_block and dirname not in b_block: c_block.append(dirname)

	for dirname in g_dirs:
		if dirname not in a_block and dirname not in b_block: c_block.append(dirname)
		
	#call check_fivematch to check the blocks
	#print 'Checking a block...'
	#check_fivematch(a_block)

	#print 'Checking b block...'
	#check_fivematch(b_block)
	
		
	print
	return a_block, b_block, c_block
            

#----------------------------------------------------------------------------------------------------#
#check the matches made in match all colours. If any colour's ccd1 if shifted more than 5 pixels
#from red, it is removed from the match list

#def check_fivematch_list(expath, match, remainder_rdirs, remainder_nbdirs, remainder_idirs, remainder_udirs, remainder_gdirs):
def check_fivematch(match):
    #check shiftx, shifty for ccd1 in each directory
    #print "Checking list to make sure shifts between ccds are <20 pixels..."
    print
    print match
    r_r_ccd = glob.glob(match[2]+'/confcorr/*_ccd1_*.fit')
    r_b_ccd = glob.glob(match[3]+'/confcorr/*_ccd1_*.fit')
    nb_ccd = glob.glob(match[5]+'/confcorr/*_ccd1_*.fit')
    i_ccd = glob.glob(match[4]+'/confcorr/*_ccd1_*.fit')
    u_ccd = glob.glob(match[0]+'/confcorr/*_ccd1_*.fit')
    g_ccd = glob.glob(match[1]+'/confcorr/*_ccd1_*.fit')

    r_r = fits.open(r_r_ccd[0])
    r_b = fits.open(r_b_ccd[0])
    nb = fits.open(nb_ccd[0])
    i = fits.open(i_ccd[0])
    u = fits.open(u_ccd[0])
    g = fits.open(g_ccd[0])

    r_r_header = r_r[0].header
    r_b_header = r_b[0].header
    nb_header = nb[0].header
    i_header = i[0].header
    u_header = u[0].header
    g_header = g[0].header
    
    r_r_wcs = wcs.WCS(r_r_header)
    r_b_wcs = wcs.WCS(r_b_header)
    nb_wcs = wcs.WCS(nb_header)           
    i_wcs = wcs.WCS(i_header)
    u_wcs = wcs.WCS(u_header)
    g_wcs = wcs.WCS(g_header)

    r_r_img = r_r[0].data
    r_b_img = r_b[0].data
    nb_img = nb[0].data 
    i_img = i[0].data   
    u_img = u[0].data
    g_img = g[0].data       

    
    cx = r_r_img.shape[0]/2
    cy = r_r_img.shape[1]/2
    cra, cdec = r_r_wcs.wcs_pix2world(cx, cy, 1) #colour1 as reference
    nbx, nby = nb_wcs.wcs_world2pix(cra, cdec, 1) #pixels in colour2 with same wcs as cra and cdec in colour1
    ix, iy = i_wcs.wcs_world2pix(cra, cdec, 1)
    ux, uy = u_wcs.wcs_world2pix(cra, cdec, 1)
    gx, gy = g_wcs.wcs_world2pix(cra, cdec, 1)
    rx, ry = r_b_wcs.wcs_world2pix(cra, cdec, 1)
    
    print
    nbshiftx = nbx-cx
    nbshifty = nby-cy
    print "nbshiftx, nbshifty: ", nbshiftx, nbshifty
    ishiftx = ix-cx
    ishifty = iy-cy
    print "ishiftx, ishifty: ", ishiftx, ishifty
    ushiftx = ux-cx
    ushifty = uy-cy
    print "ushiftx, ushifty: ", ushiftx, ushifty
    gshiftx = gx-cx
    gshifty = gy-cy
    print "gshiftx, gshifty: ", gshiftx, gshifty
    rshiftx = rx-cx
    rshifty = ry-cy
    print "rbshiftx, rbshifty: ", rshiftx, rshifty
        
    r_r.close()
    r_b.close()
    nb.close()
    i.close()
    u.close()
    g.close()
    """    
    if abs(nbshiftx)>20 or abs(nbshifty)>20:
        #shifts are too large pairing must be wrong
        print "Shift > 20 between: "
        print match[0]
        print match[1]
        print "Directories must be incorrectly paired."
       # remainder_rdirs.add(match[0])
       # remainder_nbdirs.add(match[1])
        remove = True
    if abs(ishiftx)>5 or abs(ishifty)>5: 
        print "Shift > 20 between: "
        print match[0]
        print match[2]
        print "Directories must be incorrectly paired."
      #  remainder_rdirs.add(match[0])
      #  remainder_idirs.add(match[2])
        remove = True
    if abs(ushiftx)>20 or abs(ushifty)>20:
        print "Shift > 20 between: "
        print match[0]
        print match[3]
        print "Directories must be incorrectly paired."
      #  remainder_rdirs.add(match[0])
      #  remainder_udirs.add(match[3])
        remove = True
    if abs(gshifty)>20 or abs(gshifty)>20: 
        #shifts are too large pairing must be wrong
        print "Shift > 20 between: "
        print match[0]
        print match[4]
        print "Directories must be incorrectly paired. "    
      #  remainder_rdirs.add(match[0])
      #  remainder_gdirs.add(match[4])
        remove = True
            
    if remove == True:
        print "Deleting match"
        match = []
   """

            
    print
    print "Checking completed."
    print
    return match

        
