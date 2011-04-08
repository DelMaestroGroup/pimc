# getresults.py
# Adrian Del Maestro
# 01.l5.2010
# 
# Get a bunch of results from my various clusters

import os,sys
from optparse import OptionParser

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

	# setup the command line parser options 
	parser = OptionParser() 
	parser.add_option("-T", "--temperature", dest="T", type="float", \
			help="simulation temperature in Kelvin") 
	parser.add_option("-V", "--volume", dest="V", type="float",\
			help="volume in Angstroms^d") 
	parser.add_option("-a", "--potential-width", dest="a", type="string",\
			help="width of gaussian potential") 

	(options, args) = parser.parse_args() 

	# We get the wildcard string
	if options.T and not options.V:
		wc = '*%06.3f*' % options.T
	elif options.V and not options.T:
		wc = '*%07.3f*' % options.V
	elif options.V and options.T:
		wc = '*%06.3f-%07.3f*' % (options.T,options.V)
	else:
		parser.print_help()
		sys.exit()

	# Create the rsync command 
	outDir = 'Projects/PIMC/1DBOSEGAS/LL/a_eq_' + options.a + '/OUTPUT/'
	hosts = ('orcinus.westgrid.ca:~/'+outDir,\
			'glacier.westgrid.ca:~/'+outDir,\
			'brown.sharcnet.ca:/work/agdelma/'+outDir,\
			'whale.sharcnet.ca:/work/agdelma/'+outDir)
	rsyncCmd = 'rsync -avz --exclude=*state* '

	# get the files
	for host in hosts:
		print '\nFetching results from: ', host.split('.')[0], '.....\n'
		os.system(rsyncCmd + host + wc + ' .')

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
	main()

