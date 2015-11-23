import os, glob
from pyraf import iraf

# takes a number of list files, creates folders in the working directory, and copies the files in each list to the new folder
prefixDict = {'n1':'17','n2':'18','n3':'19','n4':'20','n5':'21'}
linePrefix1 = 'data/201205'
linePrefix3 = '/reduce/'
iraf.copy.setParam('verbose','yes')
for file_ in glob.glob('*.list'):
	with open(file_) as f:
		for line in f:
			if line.startswith('#'):
				print "creating folder", line[2:-1]
				os.mkdir(line[2:-1])
				folder = line[2:-1]
			else :
				cpfile = linePrefix1+prefixDict[line[0:2]]+linePrefix3+line[0:-1]+'.fits'
				dest = folder+os.sep+line[0:-1]+'.fits'
				iraf.copy(cpfile,dest)
				
				
	