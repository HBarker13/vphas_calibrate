#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#create r-v (x-axis) vs r-Halpha (y-axis) colour-colour diagram from stars. Uses bandmerged catalgoues 
#use only the best stars ie. good photometry, magnitudes ~13 to ~19
#calculate the u bannd magnitude shift needed to put the maximum number of stars between the MS and G0V lines
#update catalogue magnitudes and write the shifts and errors to file
#r2 (blue blcok) mags are used ini the calibration


from astropy.io import fits
import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
from scipy import interpolate
import operator
import argparse
from scipy.optimize import curve_fit
import math


import make_lists
#change font size
import matplotlib
matplotlib.rcParams.update({'font.size': 20})


def append_table(table, name, arr, d_type):
    arr = np.asarray(arr)
    dtype = d_type
    newdtype = np.dtype(table.dtype.descr + [(name, d_type)])
    newtable = np.empty(table.shape, dtype=newdtype)
    for field in table.dtype.fields:
        newtable[field] = table[field]
    newtable[name] = arr
    return newtable

    
parser = argparse.ArgumentParser(description="Collect calibration inputs")
parser.add_argument('-v','--vphas_num', help="vphas pointing number", required=True)
parser.add_argument('-b','--block', help="vphas offset block", required=True)
args = parser.parse_args()

ex_path = os.getcwd() + '/vphas_' + args.vphas_num + '_ex'
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)
vphas_filternames = {'u':0, 'g':1, 'r':2, 'r2':3, 'i':4, 'NB':5}





#block_choice = raw_input('Block a or b? ')
block_choice = args.block
if block_choice=='a':
	block = a_block
	merged_fpath = os.getcwd()+'/a_block_merged_cat.fits'
elif block_choice=='b':
	block = b_block
	merged_fpath = os.getcwd()+'/b_block_merged_cat.fits'
elif block_choice=='c': 
	block = c_block
print 'Block', block_choice
print

if not os.path.exists(merged_fpath):
	print 'File does not exists:'
	print merged_fpath
	raw_input('')

openfile = fits.open(merged_fpath)
table = openfile[1].data
names = table.dtype.names
openfile.close()




"""
NEEDS REDOING
#CS colours at Rv = 3.1
#E(B-V): u-g, g-r, r-i
CS = {0:[-1.598, -0.315, -0.169], 0.2:[-1.385, -0.087, -0.048], 0.4:[-1.163, 0.132, 0.074], 0.6:[-0.935, 0.347, 0.197], 0.8:[-0.703, 0.558, 0.319], 1.0:[-0.466, 0.766, 0.442], 1.2:[-0.224, 0.970, 0.566], 1.4:[0.023, 1.170, 0.689], 1.6:[0.273, 1.367, 0.813], 1.8:[0.528, 1.561, 0.937] }

CS_u_min_g = [CS[c][0] for c in CS]
CS_g_min_r = [CS[c][1] for c in CS]
CS_r_min_i = [CS[c][2] for c in CS]
"""



#A0 synthetic colours for MS from vphas table a2. a_0 = reddening at 5500angstroms
#A_0: u-g, g-r, r-i, r-ha
#Rv=3.1
A0 = {0:[-0.053, -0.005, -0.009, -0.005], 2:[0.675, 0.780, 0.418, 0.133], 4:[1.431, 1.540, 0.833, 0.246], 6:[2.153, 2.277, 1.238, 0.334], 8:[2.514, 2.995, 1.633, 0.399], 10:[float('nan'), float('nan'), 2.020, 0.441]}

#Rv=2.5
#A0 = {0:[-0.053, -0.005, -0.009, -0.005], 2:[0.957, 0.940, 0.474, 0.155], 4:[1.987, 1.841, 0.939, 0.280], 6:[2.792, 2.706, 1.388, 0.372], 8:[2.489, 3.543, 1.823, 0.433], 10:[float('nan'), float('nan'), 2.245, 0.467]}

#Rv=3.8
#A0 = {0:[-0.053, -0.005, -0.009, -0.005], 2:[0.522, 0.652, 0.373, 0.116], 4:[1.118, 1.295, 0.749, 0.218], 6:[1.713, 1.925, 1.118, 0.300], 8:[2.199, 2.542, 0.481, 0.365], 10:[float('nan'), float('nan'), 1.839, 0.412]}


A0_r_min_i = [A0[a][2] for a in A0]
A0_g_min_r = [A0[a][1] for a in A0]
A0_u_min_g = [A0[a][0] for a in A0]
A0_r_min_ha = [A0[a][3] for a in A0]




#A2 synthetic colours for MS from vphas table a2. a_0 = reddening at 5500angstroms
#A_0: u-g, g-r, r-i, r-ha
#Rv=3.1
A2 = { 0: [0.021, 0.025, 0.006, -0.004], 2:[0.749, 0.809, 0.432, 0.134 ], 4:[1.505, 1.568, 0.874, 0.247], 5:[2.217, 2.304, 1.251, 0.334], 8:[2.538, 3.021, 1.646, 0.397], 10:[float('nan'), float('nan'), 2.033, 0.439 ] } 

A2_g_min_r = [A2[a][1] for a in A2]
A2_u_min_g = [A2[a][0] for a in A2]



#A3 synthetic colours for MS from vphas table a2. a_0 = reddening at 5500angstroms, Rv=3.1
#A_0: u-g, g-r, r-i, r-ha
A3 = {0:[0.038, 0.059, 0.021, -0.008], 2:[0.771, 0.840, 0.446, 0.130], 4:[1.531, 1.597, 0.861, 0.241], 6:[2.241, 2.332, 1.265, 0.328], 8:[2.541, 3.048, 1.660, 0.391], 10:[float('nan'), float('nan'), 2.047, 0.432]}

A3_g_min_r = [A3[a][1] for a in A3]
A3_u_min_g = [A3[a][0] for a in A3]




#G0V synthetic colours for MS with Rv=3.1 from vphas table a2
#u-g, g-r
G0V = {0:[-0.001,0.630], 2:[0.779, 1.388], 4:[1.566, 2.124], 6:[2.201,2.842], 8:[1.754, 4.081]}

#Rv=2.5
#G0V = {0:[-0.001,0.630], 2:[1.075, 1.540], 4:[2.117, 2.412], 6:[3.257, 1.754], 8:[2.140, 3.547]}

#Rv=3.8
#G0V = {0:[-0.001,0.630], 2:[0.616, 1.266], 4:[1.244, 1.890], 6:[1.822, 2.503], 8:[2.122, 3.107]}


#unreddened MS: A_0 = 0 vals from vphas appendix tables(Rv=2.5, but it doesn't matter)
#spec_type: r-i, r-nb, u-g, g-r
main_sequence = {'06V':[-0.145, 0.071, -1.494, -0.299], 
		'O8V':[-0.152, 0.055, -1.463, -0.287],
		'O9V':[-0.153, 0.049, -1.446, -0.282], 
		'B0V':[-0.150, 0.054, -1.433, -0.271], 
		'B1V':[-0.136, 0.048, -1.324, -0.240], 
		'B2V':[-0.123, 0.045, -1.209, -0.218],
		'B3V':[-0.104, 0.044, -1.053, -0.186],
		'B5V':[-0.077, 0.039, -0.828, -0.139],
		'B6V':[-0.068, 0.036, -0.728, -0.121],
		'B7V':[-0.057, 0.029, -0.580, -0.100],
		'B8V':[-0.045, 0.018, -0.388, -0.076], 
		'B9V':[-0.028, 0.006, -0.198, -0.046], 
		'A0V':[-0.009, -0.005, -0.053, -0.005], 
		'A1V':[-0.003, -0.003, -0.019, 0.005], 
		'A2V':[0.006, -0.004, 0.021, 0.025],
		'A3V':[0.021, -0.008, 0.038, 0.059], 
		'A5V':[0.051, 0.005, 0.067, 0.125], 
		'A7V':[0.083, 0.027, 0.044, 0.199], 
		'F0V':[0.149, 0.084, -0.026, 0.329], 
		'F2V':[0.177, 0.109, -0.049, 0.387], 
		'F5V':[0.225, 0.149, -0.066, 0.495], 
		'F8V':[0.259, 0.173, -0.040, 0.576], 
		'G0V':[0.280, 0.188, -0.001, 0.630], 
		'G2V':[0.295, 0.197, 0.042, 0.670], 
		'G5V':[0.327, 0.217, 0.162, 0.756], 
		'G8V':[0.358, 0.233, 0.355, 0.845], 
		'K0V':[0.385, 0.245, 0.451, 0.904], 
		'K1V':[0.399, 0.251, 0.523, 0.939],
		'K2V':[0.415, 0.258, 0.602, 0.978], 
		'K3V':[0.445, 0.270, 0.756, 1.049],
		'K4V':[0.464, 0.278, 0.841, 1.092], 
		'K5V':[0.521, 0.302, 1.064, 1.198], 
		'K7V':[0.721, 0.390, 1.364, 1.386], 
		'M0V':[0.787, 0.413, 1.348, 1.394], 
		'M1V':[0.931, 0.470, 1.311, 1.422], 
		'M2V':[1.111, 0.526, 1.238, 1.425]}  


G0V_u_min_g = [G0V[a][0] for a in G0V]
G0V_g_min_r = [G0V[a][1] for a in G0V]


ms_u_min_g = [main_sequence[sp_type][2] for sp_type in main_sequence]
ms_g_min_r = [main_sequence[sp_type][3] for sp_type in main_sequence]
    
    
    
    
#remove objs with errors flagged	
for i in range(1,6):
	table = table[table['error_bit_flag_'+str(i)]==0]
	

	
u_appendix ='1'		
g_appendix ='2'
r_appendix = '3'
r2_appendix = '4'
		



#arr to store magnitude shifts
shifts_arr = []



for ap_rad in range(2,6):
	ap_name = 'Aper_flux_'+str(ap_rad)
	print ap_name
	

	#synthetic tracks are in vega, so use vega
	u_mags = table[ap_name+'_mag_'+u_appendix]

	
	#aperture corrected mags
	if ap_name+'_corr_'+g_appendix in names:
		g_mags = table[ap_name+'_corr_'+g_appendix]	
		
	else:
		print 'No corrected g mags'
		sys.exit()		
	
	if ap_name+'_corr_'+r_appendix in names:
		r_mags = table[ap_name+'_corr_'+r_appendix]
	else:
		print 'No corrected r mags'
		sys.exit()	

	if ap_name+'_corr_'+r2_appendix in names:
		r2_mags = table[ap_name+'_corr_'+r2_appendix]
	else:
		print 'No corrected r2 mags'
		sys.exit()
		
		
		
			
	#if there all the corrected mags are NaN, there are no(t enough)
	#stars to calibrate to so use raw mags	
	apcor_fail_flag = False
	
	if all(np.isnan(val) for val in g_mags) or all(np.isnan(val) for val in r_mags) or all(np.isnan(val) for val in r2_mags) :
		g_mags = table[ap_name+'_mag_'+g_appendix]
		r_mags = table[ap_name+'_mag_'+r_appendix]
		r2_mags = table[ap_name+'_mag_'+r2_appendix]
				
		apcor_fail_flag=True
		
		
	#create a directory to store images in
	if apcor_fail_flag==True:
		img_dir = os.getcwd()+'/u_fitting_apcor_fail'
	else:
		img_dir = os.getcwd()+'/u_fitting'
	if not os.path.exists(img_dir):
		os.makedirs(img_dir)
	
		
		
		
		
		

	#remove high/low mags
	u_mags = [line if 13<line<19 else float('nan') for line in u_mags]
	g_mags = [line if 13<line<19 else float('nan') for line in g_mags]
	#r_mags = [line if 13<line<19 else float('nan') for line in r_mags]
	r2_mags = [line if 13<line<19 else float('nan') for line in r2_mags]




	#calculate colours
	g_min_r = np.subtract(g_mags, r2_mags)
	u_min_g = np.subtract(u_mags, g_mags)      

	
 
	colours = [[line[0], line[1]] for line in zip(g_min_r, u_min_g) if ~np.isnan(line[0]) and ~np.isnan(line[1])]
	g_min_r = [line[0] for line in colours]
	u_min_g = [line[1] for line in colours]



	if len(g_min_r)==0:
		print 'Length of g-r list is 0'
		raw_input('Press any key to continue')
		
	if len(u_min_g)==0:
		print 'Length of u-g list is 0'
		raw_input('Press any key to continue')


	all_hist_sum = []
	#all the u band magnitude shifts to try.
	#NB This is a large range because the count number is used to plot a gaussian, from which the optimal
	# u band shift and error are calculated
	u_shifts = np.linspace(-2.0, 1.0, 3001)


	#shift the u band magnitudes until the main body of stars lies between the G0V and MS lines
	print 'Shifting u magnitudes'
	for u_shift in u_shifts:

		#print 'Shifting by:', u_shift
		shifted_u_min_g = [val+u_shift for val in u_min_g]
		
	
		#testing plot to see what area of the histogram is being looked at
		#plt.figure()
		#plt.subplots_adjust(left=None, bottom=0.13, right=None, top=None, wspace=None, hspace=None) #stop labels being clipped
		#plt.xlim(-0.5, 3.0)
		#plt.ylim(-1.4, 3)
		#plt.xlabel('g-r')
		#plt.ylabel('u-g')
	

		#create a 2D histogram of all the stars in the pointing
		nxbins = int(max(g_min_r)/0.017)
		nybins = int(max(shifted_u_min_g)/0.025)
		Hist, xedges, yedges = np.histogram2d(g_min_r , shifted_u_min_g, bins=(nybins,nxbins))
	

		#invert
		Hist = np.rot90(Hist)
		Hist = np.flipud(Hist)
	
	
		#mask out 0 values	
		Hist = np.where(Hist>0, Hist, float('NaN')) 
		extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]


					
		#min/maximum possible colours values the main cluster of stars occupies, and how many bins on 		each axis
		xmin = min(G0V_g_min_r)
		xmax = 1.4
		cut_xedges = [val for val in xedges if xmin<val<xmax]

		ymin = -0.5
		ymax = 1.0
		cut_yedges = [val for val in yedges if ymin<val<ymax]


      
		#smooth G0V line
		x = np.array(G0V_g_min_r)
		y = np.array(G0V_u_min_g)
		G0V_func = interpolate.interp1d(x,y)
		G0V_xsmooth = cut_xedges
		G0V_ysmooth = G0V_func(G0V_xsmooth)


		#smooth ms line
		x = np.array(ms_g_min_r)
		y = np.array(ms_u_min_g)
		MS_func = interpolate.interp1d(x,y)
		MS_xsmooth = cut_xedges
		MS_ysmooth = MS_func(MS_xsmooth)
	
	
	
	
		#smooth G0V plot
		#plt.plot(G0V_xsmooth, G0V_ysmooth, 'k--',)
		#plt.annotate('G0V', xy=(2.2, 1.5))
	        
		#smooth ms plot
		#plt.plot(MS_xsmooth, MS_ysmooth, 'k-')
		#plt.annotate('MS', xy=(-0.45, 0.0))

	

		#sum the bin values of all the histogram bins between these two lines   
		hist_sum = 0

		#G0V_xsmooth and MS_xsmooth are identical, and have values taken from xedges
		for i, x_pixel in enumerate(MS_xsmooth):
	
			#G0V is above MS on the plot, but it is inverted, so G0V has more negative u-g (y-axis) values
			upper_y_pixel = MS_ysmooth[i]
			lower_y_pixel = G0V_ysmooth[i]
		

			if lower_y_pixel==upper_y_pixel: continue
	
	

			#work along each row of the histogram (the xaxis) with the same y value
			for y_ind, (y_pixel, line) in enumerate(zip(yedges, Hist)):
		
	
				#see if the y pixel is between the two lines
				if lower_y_pixel < y_pixel < upper_y_pixel:
		
					#see if the xpixel is within the length of the interpolated G0V and MS 	lines
					for x_ind, (xedge, hist_bin) in enumerate(zip(xedges, line)):
						if x_pixel==xedge:

							#check the area of the histogram being used
							#Hist[y_ind][x_ind] = 100
						
							if not np.isnan(hist_bin):	
								hist_sum+=hist_bin

							
							
		#Hist = np.flipud(Hist)						
		#plt.imshow(Hist, extent=extent, cmap='autumn')		
		#plt.gca().invert_yaxis()	
		#plt.show()					
				
					
		all_hist_sum.append( hist_sum )
							

	#choose the u_shift with the maximum hist_sum. This is the initial guess
	#of the mean of the gaussian
	index, max_sum = max(enumerate(all_hist_sum), key=operator.itemgetter(1))
	shift = u_shifts[index]
	
	
	
	
		
	#fit a Gaussian to the shift vs counts to estimate the optimal shift and error
	print 'Fitting gaussian'
	def gaus(x,a,x0,sigma):
		g = [ a  * math.exp( -( val -x0)**2 / (2*sigma**2) )	for val in x]   
		return g 		
	
	
	x = u_shifts
	y = all_hist_sum
	
	#The wings of the plot aren't very gaussian, so x axis values will have to be removed
	#if index corresponding to the u shift is to the right of the centre
	#Remove the edge values of the array so there is an equal number of points either side of the peak
	if len(u_shifts)-index < index:
		half_len = len(u_shifts)-index
	else:
		half_len = index
		

	x = [ val for i,val in enumerate(x) if i in range( index-half_len, index+half_len) ]
	y = [ val for i,val in enumerate(y) if i in range( index-half_len, index+half_len) ]
		

	#Remove values from the beginning and end of the array until the best Gaussian fit is found
	#The loop stops when the "chi-squared" between the actual points and the fit Gaussian is minimised
	#Here, chi = chi_squared / len(x) . Just using chi-squared, all the points will be removed - we 
	#still want something gaussian-ish to fit to
	prev_chi = None
	gaussian_loop = True
	gaussian_failed = False
	while gaussian_loop == True:
	
		#remove the first and last values
		x = x[1:-1]
		y = y[1:-1]
		
		if len(x)<20:
			print 'Gaussian fit did not converge'
			gaussian_failed = True
			raw_input('Press any key to continue')
			gaussian_loop = False
			continue
	
		mean = shift
		sigma = math.sqrt( sum( y*(x-mean)**2 ) / sum(y) )
	
		#p0 = [ a, mean, sigma ] where a is the Gaussian vertical scaling, and the mean is the centre of the distribution
		#already calculated to be approximately equal to the u shift that maximises the number of counts  
		popt,pcov = curve_fit(gaus,x,y,p0=[max(y), shift, sigma])
		gaus_fit = gaus(x, *popt)

		
		
		
		#calculate the "chi-squared" between the actual values and the fit to see if the loop
		#needs to continue. Once the chi squared stops increasing AND the peak of the gaussian is within 
		#5% of the maximum number of counts, continue.
	
		from scipy.stats import chisquare
		chisq = chisquare(y, gaus_fit)
		chi = chisq[0]/len(x)


		#plot the gaussian fit
		#plt.figure()
		#plt.plot(x, y, 'bo')
		#plt.plot(x, gaus_fit, 'r')
		#plt.show()
			
		
		if prev_chi == None:
			prev_chi = chi
			prev_popt = popt
			prev_pcov = pcov
			continue
			
		if chi < prev_chi:
			prev_chi = chi
			prev_popt = popt
			continue
			
		if chi > prev_chi and (0.95*max_sum)<popt[0]<(1.05*max_sum):
		#if chi > prev_chi and (0.9*max_sum)<popt[0]<(1.1*max_sum):
			gaussian_loop = False
			print 'Converged'
						
			#plot the gaussian fit
			plt.figure()
			plt.plot(x, y, 'bo')
			plt.plot(x, gaus_fit, 'r')
			plt.show()
			
			continue
			
			
	# don't adjust the magnitudes if the fitting failed, but make a plot to see what's
	# going on and add corr mags to the catalogue
	if gaussian_failed == True:
	
		x = u_shifts
		y = all_hist_sum
		
		
		plt.figure()
		plt.plot(x, y, 'bo')
		plt.xlabel('u shift')
		plt.ylabel('Locus counts')
		#plt.show()
	
		savepath = img_dir + '/'+block_choice+'_Ap'+str(ap_rad)+'_FAILED.png'
		plt.savefig(savepath)
		
		
		
				
		openfile = fits.open(merged_fpath, mode='update')
		block_tab = openfile[1].data
		print 'Updating', merged_fpath

		#new_u = [float('nan') for line in u_mags]
		#block_tab[ap_name+'_corr_'+u_appendix] = new_u
			
		if not ap_name+'_corr_'+u_appendix in block_tab.dtype.names:
			new_u = [float('nan') for line in u_mags]
			block_tab = append_table( block_tab, ap_name+'_corr_'+u_appendix, new_u, '>f4')	
	
		openfile[1].data = block_tab
		openfile.close()

		
		continue


	
	
	#popt = [a, mean, sigma]
	optimal_u_shift = prev_popt[1] #=mean of gaussian fit
	optimal_counts = prev_popt[0] #number of counts corresponding to the optimal u shift. ie. the maximum of the gaussian
	
	#Noise on the highest number of counts = sqrt(counts)
	opt_counts_noise = math.sqrt(optimal_counts)




	#interpolate the gaussian fit and find the value along it that corresponds to the gaussian peak minus the noise
	#also interpolate ushift with the same increment, and find the u shift that corresponds to the noise count
	#this will be the error on the u shift
	gaus_interp = interpolate.interp1d(x, gaus_fit)
	x_new = np.linspace( min(x), max(x), 400 )
	gaus_new = gaus_interp(x_new)
	
	
	
	#u shift corresponding to the low noise counts on each side of the peak
	lowest_counts = optimal_counts - opt_counts_noise
	
	#left side of the Gaussian peak = shits more negative than the optimal value
	left_peak = [yval for xval, yval in zip(x_new, gaus_new) if xval < optimal_u_shift]
	nearest_counts_val = min( left_peak, key=lambda c: abs( c - lowest_counts) ) #nearest count value (y axis)
	low_ushift = [ xval for xval, yval in zip(x_new, gaus_new) if yval==nearest_counts_val] #nearest ushift corresponding to this (x axis)
	low_ushift = low_ushift[0]

		
	#right side of the peak = u shift values more positive than the optimal value
	right_peak = [yval for xval, yval in zip(x_new, gaus_new) if xval > optimal_u_shift]
	nearest_counts_val = min( right_peak, key=lambda c: abs( c - lowest_counts) )
	high_ushift = [ xval for xval, yval in zip(x_new, gaus_new) if yval==nearest_counts_val]
	high_ushift = high_ushift[0]


	#average the lower and upper low shifts to get the error	
	shift_err = ((optimal_u_shift - low_ushift) + (high_ushift - optimal_u_shift) ) / 2
	
	
	#Alternative error calculation
	#For a Gaussian, the error on the mean is sigma / sqrt(N)
	#shift_err = prev_popt[2] / math.sqrt( optimal_counts )		
	
		
		
	#error on the fit of the mean from scipy
	perr = np.sqrt( np.diag(prev_pcov) )
	fitting_err = perr[1] 
		
		
		
	#There will also be some systematic errors that i can't quantifty
	#eg. Whether the lines from Drew (so the stellar model atmospheres and reddening law) are good
	# Whether the stars in the field are typical	
		
	
	tot_err = math.sqrt( fitting_err**2 + shift_err**2 )	


	#print 'Initial', shift
	#print 'Gaussian mean', optimal_u_shift
	#print 'Gaussian sigma', prev_popt[2]
	#print 'Err on Gaussian mean', tot_err



	#a plot to see what's going on and check the fit is good
	plt.figure()
	plt.plot(x, y, 'bo')
	plt.plot(x_new, gaus_new, 'r-', label='Gaussian')
	plt.xlabel('u shift')
	plt.ylabel('Locus counts')
	plt.legend(loc='best')
	#plt.show()
	
	savepath = img_dir + '/'+block_choice+'_Ap'+str(ap_rad)+'.png'
	plt.savefig(savepath)






	#make the final plot with the shifted u magnitudes to make sure it looks ok						
	plt.figure()
	plt.subplots_adjust(left=None, bottom=0.13, right=None, top=None, wspace=None, hspace=None) #stop labels being clipped

	
	best_u_min_g = [val + optimal_u_shift for val in u_min_g]
	
	#check the histogram plots the same thing as a standard plot
	#plt.plot(g_min_r, best_u_min_g, 'o')
	
	
	#create a 2D histogram of all the stars in the pointing
	nxbins = int(max(g_min_r)/0.017)
	nybins = int(max(best_u_min_g)/0.025)
	Hist, xedges, yedges = np.histogram2d(g_min_r , best_u_min_g, bins=(nybins,nxbins))


	#invert: haven't flipup because we later flip the y axis
	#DO NOT comment it out in the testing histogram, the pixels will be out of order. 
	#In that histogram, flipud is used twice
	Hist = np.rot90(Hist)
	#Hist = np.flipud(Hist)

	
	#mask out 0 values	
	Hist = np.where(Hist>0, Hist, float('NaN')) 
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

	plt.clf()
	plt.imshow(Hist, extent=extent, cmap='autumn')	
						

	plt.xlim(-0.5, 3.0)
	plt.ylim(-1.4, 3)
    
	#invery y axis
	plt.gca().invert_yaxis()

	plt.xlabel('g-r')
	plt.ylabel('u-g')

	#smooth G0V plot
	x = np.array(G0V_g_min_r)
	y = np.array(G0V_u_min_g)
	func = interpolate.interp1d(x,y)
	x_smooth = np.linspace(x.min(), x.max(), 300)
	y_smooth = func(x_smooth)
	plt.plot(x_smooth, y_smooth, 'k--',)
	plt.annotate('G0V', xy=(2.2, 1.5))
        

	#smooth ms plot
	x = np.array(ms_g_min_r)
	y = np.array(ms_u_min_g)
	func = interpolate.interp1d(x,y)
	x_smooth = np.linspace(x.min(), x.max(), 300)
	y_smooth = func(x_smooth)
	plt.plot(x_smooth, y_smooth, 'k-')
	plt.annotate('MS', xy=(-0.45, 0.0))

	

	#smooth A3 line
	x = np.array(A3_g_min_r)
	y = np.array(A3_u_min_g)
	func = interpolate.interp1d(x,y)
	xsmooth = np.linspace(x.min(), x.max(), 300)
	ysmooth = func(xsmooth)
	plt.plot(x, y, 'k--')
 	plt.annotate('A3V', xy=(2.2, 2.5))
 	

 	
 	#save
	savepath = img_dir +'/'+block_choice+'_'+ap_name+'.png'
	plt.savefig(savepath)
	print 'Saved', savepath
	#plt.show()



	#keep track of the calculated shifts and errors
	shifts_arr.append( [str(ap_rad), str(optimal_u_shift), str(shift_err) ] )






	#add corrected u_mags to fits file
	openfile = fits.open(merged_fpath, mode='update')
	block_tab = openfile[1].data
	print 'Updating', merged_fpath
	
	if not ap_name+'_corr_'+u_appendix in block_tab.dtype.names:
		u_mags = block_tab[ap_name+'_mag_'+u_appendix]
		new_u = [line + optimal_u_shift for line in u_mags]
		block_tab = append_table( block_tab, ap_name+'_corr_'+u_appendix, new_u, '>f4')
	
	else:
		u_mags = block_tab[ap_name+'_mag_'+u_appendix]
		new_u = [line + optimal_u_shift for line in u_mags]
		block_tab[ap_name+'_corr_'+u_appendix] = new_u
	
	openfile[1].data = block_tab
	openfile.close()


	
	catpath = glob.glob(block[0]+'/catalogues/*_cat.fits',)
	openfile = fits.open(catpath[0], mode='update')
	print 'Updating', catpath[0]

	for ccdnum in range(1,33):
		#print 'ccd', str(ccdnum)
		ccd_tab = openfile[ccdnum].data
		
		if not ap_name+'_corr' in ccd_tab.dtype.names:
			u_mags = ccd_tab[ap_name+'_mag']
			new_u = [line + optimal_u_shift for line in u_mags]
			ccd_tab = append_table(ccd_tab, ap_name+'_corr', new_u, '>f4')

		else:
			u_mags = ccd_tab[ap_name+'_mag']
			new_u = [line + optimal_u_shift for line in u_mags]
			ccd_tab[ap_name+'_corr']=new_u

		openfile[ccdnum].data = ccd_tab
	openfile.close()
	print







#write the shifts to a file
savepath = os.getcwd()+'/u_corrections_'+block_choice+'.txt'
with open(savepath, 'w+') as f:
	for line in shifts_arr:
		f.write(line[0]+'	'+line[1]+'	'+str(line[2])+'	\n')
print 'Saved', savepath
print	
















