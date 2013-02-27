"""Plot Tensor Estimator

Description:
  Performs a density plot of a tensor estimator.  Currently defaults to 3d
  
Usage:
  plot_tensor_estimator.py [--estimator=<estimator> --Lx=<Lx> --Ly=<Ly> --Lz=<Lz>] <input-file>

  plot_tensor_estimator.py -h | --help 

Options:
  -e <estimator>, --estimator=<estimator> The name of an estimator to be plotted
  --Lx=<Lx>                               The spatial extent in the x-direction
  --Ly=<Ly>                               The spatial extent in the y-direction
  --Lz=<Lz>                               The spatial extent in the z-direction
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
import loadgmt,kevent
import itertools

#cmap = loadgmt.getcmap('dca','alarm.p1.0.1')
#cmap = loadgmt.getcmap('cb/div','Spectral_11')
#cmap = loadgmt.getcmap('cb/div','RdGy_11')
cmap = loadgmt.getcmap('gist','earth')
cmap = loadgmt.getcmap('ncl','precip3_16lev')
cmap = loadgmt.getcmap('ncl','temp_diff_18lev')

# ===================================================================
def main():

    # Read in the command line arguments
    args = docopt(__doc__)
    fileName = args['<input-file>']
    estName = args['--estimator']

    # Open up the tensor file, and determine the number of grid boxes in each
    # dimension and the column headers
    with open(fileName,'r') as inFile:
        N = int(inFile.readline().split()[1])

    # Get the spatial dimensions, otherwise just use N
    L = [N,N,N]
    for n,l in enumerate(['--Lx','--Ly','--Lz']):
        if args[l]:
            L[n] = float(args[l])

    gridSize = []
    for cL in L:
        gridSize.append(cL/(1.0*N))
    dV = pl.product(gridSize)

    # Get the column headers from the data file
    headers = pimchelp.getHeadersDict(fileName,skipLines=1)

    # determine which column we will plot
    plotColumn = headers.get(estName,0)

    # Assuming a 3d data set, load the data
    data = pl.loadtxt(fileName,ndmin=2)[:,plotColumn].reshape([N,N,N])

    # assuming translational symmetry in the z-axis
    rhoxy = pl.average(data,axis=2)

    # For the particular case of A^2 we need to 

    num_rad_sep = int(N/2)
    dr = 0.5*L[0]/num_rad_sep
    rho = pl.zeros(num_rad_sep)
    counts = pl.zeros_like(rho)
    l = range(N)
    for i,j in itertools.product(l,l):
        x = -0.5*L[0] + (i+0.5)*gridSize[0]
        y = -0.5*L[1] + (j+0.5)*gridSize[1]
        r = pl.sqrt(x*x + y*y)
        n = int(r/dr)
        if n < num_rad_sep:
            if abs(rho[n]) < 1000.0: 
                if 'A^2' in estName:
                    rho[n] += rhoxy[i,j]/(r*r)
                else:
                    rho[n] += rhoxy[i,j]
                counts[n] += 1
        if 'A^2' in estName and r > 0.8:
            rhoxy[i,j] /= r*r
        elif 'A^2' in estName and r < 0.8:
            rhoxy[i,j] = 0.0
        elif 'A:rho_s' in estName and r < 1.0:
            rhoxy[i,j] = 0.0
        
    for n in range(num_rad_sep):
        if counts[n] > 0:
            rho[n] /= counts[n]

    # plot histograms in all three projections
    pl.figure(1)
    pl.imshow(rhoxy, cmap=cmap, 
              extent=[-0.5*L[0],0.5*L[0],-0.5*L[1],0.5*L[1]], 
             interpolation='None')
    pl.xlabel(r'$y\  [\AA]$')
    pl.ylabel(r'$x\  [\AA]$')
    pl.title('Particle Density Projection (X-Y)')
    pl.colorbar(shrink=0.4)

    r = pl.linspace(0,0.5*L[0],num_rad_sep)
#    ar = pl.loadtxt('r.dat')
    pl.figure(2)
    pl.plot(r,rho,'ro-', markersize=8)
#    pl.plot(ar[:,0],ar[:,1],'r-', linewidth=2)

    for i,cr in enumerate(r):
        print '%16.8E%16.8E' % (cr,rho[i])
    
    # Now we locally average the results
    loc_ave_data = pl.zeros_like(rhoxy)
    for i,j in itertools.product(l,l):
        x = -0.5*L[0] + i*gridSize[0]
        y = -0.5*L[1] + j*gridSize[1]
        r = pl.sqrt(x*x + y*y)
        if r < 0.5*L[0]:
            count = 0
            for k,l in itertools.product(range(-1,2),range(-1,2)):
                if i+k < N and j+l < N:
                    loc_ave_data[i,j] += rhoxy[i+k,j+l]
                    count += 1
            if count > 0:
                loc_ave_data[i,j] /= count
        else:
            loc_ave_data[i,j] = rhoxy[i,j]

    pl.figure(3)
    pl.imshow(loc_ave_data, cmap=cmap, 
              extent=[-0.5*L[0],0.5*L[0],-0.5*L[1],0.5*L[1]])
    pl.xlabel(r'$y\  [\AA]$')
    pl.ylabel(r'$x\  [\AA]$')
    pl.title('Particle Density Projection (X-Y)')
    pl.colorbar(shrink=0.4)
#    pl.figure(2)
#    pl.imshow(I, cmap=cmap, extent=[-0.5*L[0],0.5*L[0],-0.5*L[1],0.5*L[1]])
#    pl.xlabel(r'$y\  [\AA]$')
#    pl.ylabel(r'$x\  [\AA]$')
#    pl.title('Particle Density Projection (X-Y)')
#    pl.colorbar(shrink=0.4)
    
#    pl.figure(3)
#   # pl.imshow(R2, cmap=cmap, 
#   #           extent=[-0.5*L[0],0.5*L[0],-0.5*L[1],0.5*L[1]])
#    pl.imshow(pl.sum(data,axis=1), cmap=cmap, 
#              extent=[-0.5*L[2],0.5*L[2],-0.5*L[0],0.5*L[0]])
#    pl.xlabel(r'$z\  [\AA]$')
#    pl.ylabel(r'$x\  [\AA]$')
#    pl.title('Particle Density Projection (X-Z)')
#    pl.colorbar(shrink=0.4)
#  
#    pl.figure(4)
#    pl.imshow(pl.sum(data,axis=0), cmap=cmap,
#              extent=[-0.5*L[2],0.5*L[2],-0.5*L[1],0.5*L[1]])
#    pl.xlabel(r'$z\  [\AA]$')
#    pl.ylabel(r'$y\  [\AA]$')
#    pl.title('Particle Density Projection (Y-Z)')
#    pl.colorbar(shrink=0.4)
#    pl.savefig('plot.svg')

   
    pl.show()
          
# ===================================================================
if __name__=="__main__":
    main()
