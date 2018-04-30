import numpy as np
import os
from magfilter import magfilter

path = os.getcwd()
steps = path.split('/')
title_string = steps[-1].upper()

filter_file = os.path.dirname(os.path.abspath(__file__))+'/filter.txt'
filter_string = 'old'

dm2, mpc2, sep2, n2, s2, pct2, pct_hi2 = np.loadtxt('search_2.0.txt', usecols=(0,1,4,5,6,7,8), unpack=True)
dm3, mpc3, sep3, n3, s3, pct3, pct_hi3 = np.loadtxt('search_3.0.txt', usecols=(0,1,4,5,6,7,8), unpack=True)

ra2, dec2 = np.loadtxt('search_2.0.txt', usecols=(2,3), dtype=str, unpack=True)
ra3, dec3 = np.loadtxt('search_3.0.txt', usecols=(2,3), dtype=str, unpack=True)

close2 = np.where(sep2 < 480.)
close3 = np.where(sep3 < 480.)

maxloc2 = np.argmax(pct2)
maxloc3 = np.argmax(pct3)


print len(np.where(pct2>90)[0]), 'significant overdensities'
print len(np.where(sep2<480)[0]), 'close overdensities'
print len(np.where((pct2>90) & (sep2<480))[0]), 'significant and close overdensities'

print '======== best ==========='
dm2b, mpc2b, ra2b, dec2b, sep2b, n2b, s2b, pct2b, hi2b = magfilter(2.0, '2_0', dm2[maxloc2], '{:5.2f}'.format(dm2[maxloc2]), filter_file, filter_string)
try:
	maxclose2 = np.argmax(pct2[close2])
	print '======== close ==========='
	dm2c, mpc2c, ra2c, dec2c, sep2c, n2c, s2c, pct2c, hi2c = magfilter(2.0, '2_0', (dm2[close2])[maxclose2], '{:5.2f}'.format((dm2[close2])[maxclose2]), filter_file, filter_string)
except:
	print 'no close overdensities'

print len(np.where(pct3>90)[0]), 'significant overdensities'
print len(np.where(sep3<480)[0]), 'close overdensities'
print len(np.where((pct3>90) & (sep3<480))[0]), 'significant and close overdensities'

print '======== best ==========='
dm3b, mpc3b, ra3b, dec3b, sep3b, n3b, s3b, pct3b, hi3b = magfilter(3.0, '3_0', dm3[maxloc3], '{:5.2f}'.format(dm3[maxloc3]), filter_file, filter_string)
try:
	maxclose3 = np.argmax(pct3[close3])
	print '======== close ==========='
	dm3c, mpc3c, ra3c, dec3c, sep3c, n3c, s3c, pct3c, hi3c = magfilter(3.0, '3_0', (dm3[close3])[maxclose3], '{:5.2f}'.format((dm3[close3])[maxclose3]), filter_file, filter_string)

except:
	print 'no close overdensities'

print sep3b-sep2b, s3b-s2b, pct3b-pct2b
print sep3c-sep2c, s3c-s2c, pct3c-pct2c