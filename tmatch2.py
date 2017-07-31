#!/usr/local/anaconda/bin/python

"""Create fits file of matching stars in apass and the vphas point source catalogue."""

import os
import sys
import stilts

#get table filepaths. sys.argv[0] = script name
t1name = sys.argv[1]
t2name = sys.argv[2]

intype1 = sys.argv[3]
radtodeg1 = sys.argv[4]

intype2 = sys.argv[5]
radtodeg2 = sys.argv[6]

ra1 = sys.argv[7]
dec1 = sys.argv[8]
ra2 = sys.argv[9]
dec2 = sys.argv[10]

match_radius = sys.argv[11]
jointype = sys.argv[12]
findtype = sys.argv[13]
savename = sys.argv[14]

print 'Creating: ', savename


t1 = stilts.tread(t1name, intype1)
t2 = stilts.tread(t2name, intype2)


if radtodeg1=='True':
	ra1 = 'radiansToDegrees('+ra1+')'
	dec1 = 'radiansToDegrees('+dec1+')'
	
if radtodeg2=='True':
	ra2 = 'radiansToDegrees('+ra2+')'
	dec2 = 'radiansToDegrees('+dec2+')'



#match to within match_radius. Keep only rows with a value from table 1 and table 2. find=a'all' keepa all matches between tables 1 and 2 (don't choose the best one), 'best' chooses the best
matched = stilts.tskymatch2(in1=t1, in2=t2, ra1=ra1, dec1=dec1, ra2=ra2, dec2=dec2, error=match_radius, join=jointype, find=findtype)
matched.write(savename)


