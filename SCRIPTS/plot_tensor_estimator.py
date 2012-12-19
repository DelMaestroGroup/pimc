"""Plot Tensor Estimator

Description:
  Performs a density plot of a tensor estimator.  Currently defaults to 3d
  
Usage:
  plot_tensor_estimator.py [--estimator=<estimator>] <input-file>

  plot_tensor_estimator.py -h | --help 

Options:
  -e <estimator>, --estimator=<estimator> The name of an estimator to be plotted
  -h --help                               Show this screen.
"""

# ===================================================================
# Plot the density of particles in a simulation cell.
#
# Author:       Max Graves
#               Adrian Del Maestro
# Date:         8-NOV-2012
# ===================================================================

import pylab as pl
import argparse
import pimchelp
from docopt import docopt

# ===================================================================
def main():

    # Read in the command line arguments
    args = docopt(__doc__)
    fileName = args['<input-file>']
    estName = args['--estimator']

    # Open up the tensor file, and determine the number of grid boxes in each
    # dimension and the column headers
    inFile = open(fileName,'r')
    N = int(inFile.readline().split()[1])
    inFile.close()

    # Get the column headers from the data file
    headers = pimchelp.getHeadersDict(fileName,skipLines=1)

    # determine which column we will plot
    if estName and estName in headers:
        col = headers[estName]
    else:
        col = 0

    # Assuming a 3d data set, load the data
    data = pl.loadtxt(fileName,ndmin=2)[:,col].reshape([N,N,N])
        
    # plot histograms in all three projections
    pl.figure(1)
    pl.imshow(pl.sum(data,axis=2))
    pl.xlabel(r'$y\  [\AA]$')
    pl.ylabel(r'$x\  [\AA]$')
    pl.title('Particle Density Projection (X-Y)')
    pl.colorbar(shrink=0.4)
    
    pl.figure(2)
    pl.imshow(pl.sum(data,axis=1))
    pl.xlabel(r'$z\  [\AA]$')
    pl.ylabel(r'$x\  [\AA]$')
    pl.title('Particle Density Projection (X-Z)')
    pl.colorbar(shrink=0.4)
   
    pl.figure(3)
    pl.imshow(pl.sum(data,axis=0))
    pl.xlabel(r'$z\  [\AA]$')
    pl.ylabel(r'$y\  [\AA]$')
    pl.title('Particle Density Projection (Y-Z)')
    pl.colorbar(shrink=0.4)

   
    pl.show()
          
# ===================================================================
if __name__=="__main__":
    main()
