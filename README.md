#----------------To calibrate a VPHAS+ pointing-----------------#


Assumes full_mosaic.py has already been done, so the files have already been downloaded and sorted_vphas.py has been run.

NB. THIS ONLY CALIBRATES APERTURES 3,4,5,6,7. THIS DOES NOT CALIBRATE THE HALPHA DATA


1. Download VPHAS+ files: download_vphas.sh

2. Sort the files using sort_vphas.py

3. Decompress all the files: decompress_vphas.py

4. Extract the files: extract_ccds.py


5. calc_ap_offset.py.   Download apass data from: https://www.aavso.org/apass  and save as a .csv file.  I use the coordinates of the vphas pointing and a 2 degree radius.  This script loops throught the apass data and removes mutliple entries of the same object. This takes a while but not doing this introduces a lot of noise in the calibration.
	
	Apass and vphas objects are then matched, with any poot photometry thrown away.  The stars are matched by matching all vphas stars within 5arcsec of an apass star and choosing the brightnest as the true match. (If there was a brighter star neaby, that one would have been chosen as the apass star). Again, this reduces noise in the calibration.  
	
	The magnitudes of the matched apass and vphas stars are calculated (in AB) for apertures 3,4,5,6 and plotted, and the average magnitude offset calculated. Creates txt files with aperture corrections and errors.
	
	NB. THE USER MUST CHECK the plots created to make sure there's enough apass stars for a good calibration. However, the error comes from the standard deviation of distance from the line, so poor fits will have a larger error.
	

6. process_catalogues.py  (also known as combine_cats.py).  Loops through all the catalogue files, adding useful header information as columns to make using the bandmerged catalogue easier. Uses the txt files of aperture corrections to correct the vpahs magntiudes. An error will be flagged for some of the apertures 

7. calc_errs.py  Calculate the errors from the photon count, zeropoint and aperture correction error, and adds columns to the catalogue files

8. bandmerge_catalogue.py  Bandmerges the catalogues for each block using Topcat. ie. Block A contains u,g,r,r2,i 

9. u_correct.py   Uses the bandmerged catalogues. Calibrate the u band magnitudes by plotting the calibrated u-g, g-r colours of the vphas stars. The user must choose the u band correction to apply so that largest density of stars falls between the MS and G0V reddening line. The corrected u band magnitudes and errors are added to the bandmerged catalogue and the single band catalogue. 


