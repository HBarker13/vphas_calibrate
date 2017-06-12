#----------------To calibrate a VPHAS+ pointing-----------------#



NB. This only calibrates apertures 3,4,5,6,7, and does not calibrate Halpha data. It assumes a complete data set has been downloaded from CASU. ie. There are 2 u, r, r2, i and 3 g and NB exposures. 
All the below scripts can be called using calibrate_wrapper.py. If you do this, be sure to check the aperture correction plots  (in /aperture_to_apass) to check there was enough vphas - apass matches for a good fit. If the fit is too poor, try commenting out some of the apass star removal criteria (use with caution!).

The scripts assume they are been run from the directory outside vphas_xxx. This can be changed by changing the 'working_path' in the scripts. They can be called using:
	python script_name.py -v xxxx
This was added incase there's more than one vphas pointing in the directory.




1. Download VPHAS+ files from CASU. This set of scripts they are put in a directory vphas_xxxx  where xxx if the four number vphas pointing number. eg. vphas_0175. Move the filelist to /vphas_xxx/file_list_xxx 

2. sort_vphas.py  Make a directory vphas_xxx_ex and sorts the raw files. This makes it easier to see which directory contains ovbservations from which band. 

3. decompress_vphas.py : Decompress all the files using imcopy

4. extract_ccds.py : Extract each of the 32 CCDs into its own file

5. reduce_apass.py   Download apass data from: https://www.aavso.org/apass  and save as a .csv file.  I use the coordinates of the vphas pointing and a 2 degree radius. Remove very faint / very bright stars in the apass catalalogue, combine observations of the same star that aren't listed with the same coordinates, and remove objects with only one observation. Creates apass_reduced.csv.


6. calc_ap_offset.py.   
	Apass and vphas objects are then matched, with any poot photometry thrown away.  The stars are matched by matching all vphas stars within 5arcsec of an apass star and choosing the brightnest as the true match. (If there was a brighter star neaby, that one would have been chosen as the apass star). This reduces noise in the calibration.  
	
	The magnitudes of the matched apass and vphas stars are calculated (in AB) for apertures 3,4,5,6 and plotted, and the average magnitude offset calculated. Creates txt files with aperture corrections and errors.
	
	NB. THE USER MUST CHECK the plots created to make sure there's enough apass stars for a good calibration. However, the error comes from the standard deviation of distance from the line, so poor fits will have a larger error.
	

7. process_catalogues.py  (previously combine_cats.py).  Loops through all the catalogue files, adding useful header information as columns to make using the bandmerged catalogue easier. Uses the txt files of aperture corrections to correct the vpahs magntiudes. An error will be flagged for some of the apertures 

Column names: ***_mag_AB = AB magnitudes using nightzpt (from the header, in AB):   -2.5*log10( counts / exptime ) + nightzpt
              ***_mag_ = vega magnitude. AB to vega conversions calculated using zeropoints from SVO. Janet Drew says the u badn ones are likely wrong by ~0.3 mag. This should be corrected for by the calibration process
              ***_corr = vega magnitudes, including the aperture correction
              ***_corr_AB = AB magnitudes, including the aperture correction:  -2.5*log10( counts / exptime ) + nightzpt - ap_correction


8. bandmerge_catalogue.py  Bandmerges the catalogues for each block using Topcat. ie. Block A contains u,g,r,r2,i 

9. u_correct_auto.py  Uses the bandmerged catalogues. Calibrates the u band by plotting (u-g) vs (g-r) plots, overlaying stellar reddening lines from Drew 2014, and shifting u band magnitudes until the most stars are between the MS and G0V lines.
The error on the u band shift is the shift within which the star count differs by less than 5%.
Results are written to file and calibrated u band magnitudes added to catalogues.

(OLD: u_correct.py   Uses the bandmerged catalogues. Calibrate the u band magnitudes by plotting the calibrated u-g, g-r colours of the vphas stars. The user must choose the u band correction to apply so that largest density of stars falls between the MS and G0V reddening line. The corrected u band magnitudes and errors are added to the bandmerged catalogue and the single band catalogue. )


10. Run calc_errs.py to calculate the total error (photon count and aperture correction errors). Created columns containing the upper and lower magnitude limits on a measurment.


11. Run bandmerge again to get catalogues with errors
