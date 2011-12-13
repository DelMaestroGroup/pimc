#!/usr/bin/python
# pimcave.py
# Adrian Del Maestro
# 07.20.2009
# 
# Reduce and average results for a single PIMC estimator data file
# supplied as an input

import os,sys
import pyutils
from optparse import OptionParser

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = OptionParser() 
    parser.add_option("-s", "--skip", dest="skip", type="int",\
            help="how many input lines should we skip?")
    parser.set_defaults(skip=0)

    # parse the command line options and get the file name
    (options, args) = parser.parse_args() 
    if len(args) < 1: 
        parser.error("need a file name")
    
    fileNames = args

    for fileName in fileNames:
        normalize = False;

        # We check to see if we are dealing with the one body density matrix
        if fileName.find('obdm') != -1:
            normalize = True

        # We count the number of lines in the estimator file to make sure we have
        # some data and grab the headers
        estFile = open(fileName,'r');
        estLines = estFile.readlines();
        numLines = len(estLines) - 2    # We expect two comment lines
        pimcid = estLines[0]
        headers = estLines[1].split()
        estFile.close()

        # If we have data, compute averages and error
        if numLines-options.skip > 0:
            estData = pyutils.loadFile(fileName)

            # Now we skip data rows to test for convergence
            for n in range(options.skip):
                estData.pop(0)

            estAve = pyutils.average(estData,1)
            estErr = pyutils.error(estData,1)
#           estErr = pyutils.bootstrap(estData,1)
            numData = len(estData)

            print pimcid, '# Number Samples %6d' %  numData
            if not normalize:
                for n,ave in enumerate(estAve):
                    if len(headers) - 1 ==  len(estAve):
                        label = headers[n+1]
                    else:
                        label = 'Col #%02d:' % n
                    print '%-16s%12.5f\t%12.5f' % (label,estAve[n],estErr[n])
            else:
                for n,ave in enumerate(estAve):
                    normAve = estAve[n]/estAve[0]
                    if abs(estAve[n]) > 1E-10:
                        normErr = (estErr[n] / estAve[n]) * normAve
                    else: 
                        normErr = 0.0;

                    if len(headers) - 1 ==  len(estAve):
                        label = headers[n+1]
                    else:
                        label = 'Col #%02d:' % n
                    print '%-16s%12.5f\t%12.5f' % (label,normAve,normErr)
    
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
