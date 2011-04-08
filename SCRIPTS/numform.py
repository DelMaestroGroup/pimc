''' Update the number format of output files. '''
import os

outType  = ['estimator','number','obdm','pair','pcycle','super','worm','log','state']

for type in outType:
	lsCommand = 'ls -1 gce-%s-*' % type
	fileNames = os.popen(lsCommand).read().split('\n')
	fileNames.pop()

	for fname in fileNames:
		splitName = fname.split('-')
		newName = 'gce-%s-%06.3f-%07.3f-%06.3f-%04d-%s' % (splitName[1],float(splitName[2]),\
				float(splitName[3]),float(splitName[4]),int(splitName[5]),splitName[6])
		os.popen('mv %s %s' % (fname,newName))
