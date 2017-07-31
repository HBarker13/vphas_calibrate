#!/usr/local/anaconda/bin/python

"""Create fits file of matching stars in apass and the vphas point source catalogue. Merges all 6 bands (two r) at once"""

import os
import sys
import stilts

sys.setrecursionlimit(2000)

#get table filepaths. sys.argv[0] = script name
t1name = sys.argv[1]
t2name = sys.argv[2]
t3name = sys.argv[3]
t4name = sys.argv[4]
t5name = sys.argv[5]
t6name = sys.argv[6]

reftab = sys.argv[7]
match_radius = sys.argv[8] #in arcsec
savename = sys.argv[9]


t1 = stilts.tread(t1name)
t2 = stilts.tread(t2name)
t3 = stilts.tread(t3name)
t4 = stilts.tread(t4name)
t5 = stilts.tread(t5name)
t6 = stilts.tread(t6name)


#match to within 0.5arcsec
#pair matching to first input table works while multimode=group causes recursion limit error
#print 'Pair matching'
#matched = stilts.tmatchn(multimode='pair', nin=6, matcher='sky', params=match_radius, in1=t1, in2=t2, in3=t3, in4=t4, in5=t5, in6=t6, values1='RA DEC', values2='RA DEC', values3='RA DEC', values4='RA DEC', values5='RA DEC', values6='RA DEC', iref=reftab)

print 'Group matching'
matched = stilts.tmatchn(multimode='group', nin=6, matcher='sky', params=match_radius, in1=t1, in2=t2, in3=t3, in4=t4, in5=t5, in6=t6, values1='radiansToDegrees(RA) radiansToDegrees(DEC)', values2='radiansToDegrees(RA) radiansToDegrees(DEC)', values3='radiansToDegrees(RA) radiansToDegrees(DEC)', values4='radiansToDegrees(RA) radiansToDegrees(DEC)', values5='radiansToDegrees(RA) radiansToDegrees(DEC)', values6='radiansToDegrees(RA) radiansToDegrees(DEC)')


print savename
matched.write(savename)


