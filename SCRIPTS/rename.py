#!/usr/bin/python
# rename.py
# Adrian Del Maestro
# 09.28.2009
# 
# Rename a series of PIMC output files based on an input file

import os,sys,glob
from optparse import OptionParser

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

	# setup the command line parser options 
	parser = OptionParser() 

	# parse the command line options and get the file name
	(options, args) = parser.parse_args() 
	if len(args) != 1: 
		parser.error("need a file name")
	
	fileName = args[0]

	# The output file types
	fileType = ['estimator','log','obdm','pair','pcycle','state','super','worm','number']

	# We parse the input file name
	fileParts = fileName.partition('estimator')

	# Now break up the data name into descriptive parts
	dataName = fileParts[2]
	dataName = dataName.rstrip('.dat')
	dataParts = dataName.split('-')

	# Get the ID number
	oldID = int(dataParts[-1])
	newID = oldID

	# Now we keep incrementing the ID number until we are sure it is unique
	while len(glob.glob('*estimator*-%09d*' % newID)) > 0:
		newID += 1

	# Create the new data name
	dataName = ''
	for lab in dataParts[:-1]:
		dataName += lab + '-'
	dataName += '%09d.dat' % newID

	for type in fileType:
		oldName = fileParts[0] + type + fileParts[2]
		newName = fileParts[0] + type + dataName
		os.popen('mv %s %s' % (oldName,newName))

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
	main()
