"""Plot Tensor Estimator

Description:
  Performs a density plot of a tensor estimator.  Currently defaults to 3d
  
Usage:
  plot_tensor_estimator.py <input-file>

  plot_tensor_estimator.py -h | --help 

Options:
  -h --help                                 Show this screen.
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
from docopt import docopt

# ===================================================================
def RebuildTensors(Nx,Ny,Nz,data):
    """
    Takes in a row-major rank-3 tensor in single array format
    and returns the rebuilt rank-3 tensor (as 3-D numpy array)
    """
    rebXY = pl.zeros([Nz,Nx,Ny])
    rebXZ = pl.zeros([Nx,Ny,Nz])
    rebYZ = pl.zeros([Ny,Nx,Nz])

    count = 0
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                rebXY[k][i][j] = data[count]
                rebXZ[i][j][k] = data[count]
                rebYZ[j][i][k] = data[count]
                count += 1

    return rebXY, rebXZ, rebYZ

# ===================================================================
def main():

    # Read in the command line arguments
    args = docopt(__doc__)
    fileName = args['<input-file>']

    # Open up the tensor file, and determine the number of grid boxes in each
    # dimension
    inFile = open(fileName,'r')
    N = int(inFile.readline().split()[1])
    inFile.close()

    # Assuming a 3d data set, load the data
    data = pl.loadtxt(fileName).reshape([N,N,N])
        
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
