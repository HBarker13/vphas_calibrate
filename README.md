#----------------To calibrate a VPHAS+ pointing-----------------#


IMPORTANTS POINTS BEFORE YOU START

- This only calibrates apertures 3,4,5,6,7. It assumes a complete data set has been downloaded from CASU. ie. There are 2 u, r, r2, i and 3 g and NB exposures. 

- The scripts assume they are been run from the directory outside vphas_xxxx, where xxxx if the vphas pointing number. This can be changed by changing the 'working_path' in the scripts. All scripts can be called using:
	python script_name.py -v xxxx
This was added incase there's more than one vphas pointing in the directory.

- In calc_ap_offset and bandmerge_catalogue.py, jystilts is used. You will need to go into the files and in the first few rows, change the jystilts and tmatch paths. The full path must be used. Jystilts can be downloaded here: 
http://www.star.bris.ac.uk/~mbt/stilts/#install
I use the standalone Jystilts.jar. 

- The calibration is based on APASS. Data must be downlaoded in advance from: https://www.aavso.org/apass  and saved as a .csv file. reduce_apass.py will look for any file in the working directory (default /vphas_xxxx ) with "*apass*.csv" in the name. I search using the vphas pointing co-ordinates and a 1.5 degree radius.

- calibrate_halpha.py is relatively new, so can only be run seperately once the bandmerged catalogues have been created. 


IF EVERYTHING IS SET UP CORRECTLY....

...All the scripts can be called using calibrate_wrapper.py. This will call them in the correct order and pass important flags. HOWEVER, if you do this, be sure to check the aperture correction plots  (in /aperture_to_apass) to check there was enough vphas - apass matches for a good fit, and keep an eye on the progress - I wrote these scripts for me, so they are not robust and do not stop if the previous one broke. The request number is the number returned by CASU when you request the data eg. 001234 (remember the leading zeros!).




IF YOU WANT TO GO ONE STEP AT A TIME...

1. download_vphas.sh : Download VPHAS+ files from CASU. You need to go into the file beforehand and change USERNAME and PASSWORD - no quotation marks are needed. The data is put into a directory named vphas_xxxx  where xxx if the four number vphas pointing number. eg. vphas_6789. To call this script on it's own, you need to pass the request number and vphas pointing number:
	> download_vphas.sh req_num pointing_num
eg. download_vphas.sh 001234 6789


2. sort_vphas.py  This makes a directory, vphas_xxx_ex, and sorts the raw files to make it easier to see which directory contains ovbservations from which band. 
Called as: sort_vphas.py -v pointing_num
eg. sort_vphas.py -v 6789


3. decompress_vphas.py : Decompress all the files using imcopy, as vphas data is Rice compressed. This step may not be necessary?

4. extract_ccds.py : Each vphas image (single) file has a CCD's data in each of the FITS extentions. This script extracts each of the 32 CCD images into their own files

5. reduce_apass.py  Removes apass datapoints with large errors, only one observation, and merges repeated measurments of the same object.


6. calc_ap_offset.py.   
	Apass and vphas objects are then matched, with any poor photometry thrown away.  The stars are matched by matching all vphas stars within 5arcsec of an apass star and choosing the brightnest as the true match. (If there was a brighter star neaby, that one would have been chosen as the apass star). This reduces noise in the calibration.  
	
	The magnitudes of the matched apass and vphas stars are calculated (in AB) for apertures 3,4,5,6 and plotted, and the average magnitude offset calculated. Creates txt files with aperture corrections and errors.
	
	NB. THE USER MUST CHECK the plots created to make sure there's enough apass stars for a good calibration. However, the error comes from the standard deviation of distance from the line, so poor fits will have a larger error.
	

7. process_catalogues.py  (previously combine_cats.py).  Loops through all the catalogue files, adding useful header information as columns to make using the bandmerged catalogue easier to use. Uses the txt files of aperture corrections to correct vphas magntiudes.  

Column names: ***_mag_AB = AB magnitudes using apasszpt (from the header, in AB):   -2.5*log10( counts / exptime ) + apasszpt
              ***_mag_ = vega magnitude. AB to vega conversions calculated using zeropoints from SVO. Janet Drew says the u band ones are likely wrong by ~0.3 mag. This should be corrected for by the u band calibration process
              ***_corr = vega magnitudes, including the aperture correction
              ***_corr_AB = AB magnitudes, including the aperture correction:  -2.5*log10( counts / exptime ) + apasszpt - ap_correction


8. bandmerge_catalogue.py  Bandmerges the catalogues for each block using Topcat. ie. Block A contains u,g,r,r2,i 

9. u_correct_auto.py  Uses the bandmerged catalogues. Calibrates the u band by plotting (u-g) vs (g-r) plots, overlaying stellar reddening lines from Drew 2014, and shifting u band magnitudes until the most stars are between the MS and G0V lines.
The error on the u band shift is the shift within which the star count differs by less than 5%.
Results are written to file and calibrated u band magnitudes added to catalogues.

10. Run calc_errs.py to calculate the total error (photon count and aperture correction errors). Created columns containing the upper and lower magnitude limits on a measurment.


11. Run bandmerge again to get catalogues with errors
