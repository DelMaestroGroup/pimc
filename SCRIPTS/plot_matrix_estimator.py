"""Plot Matrix Estimator

Description:
  Performs a density plot of a Matrix estimator.
  
Usage:
  plot_matrix.py [--Lx=<Lx> --Ly=<Ly> --rskip=<A>] <input-file>

  plot_matrix.py -h | --help 

Options:
  -e <estimator>, --estimator=<estimator> The name of an estimator to be plotted
  --rskip=<A>, -r <A>                     The radius of the core to skip
  --Lx=<Lx>                               The spatial extent in the x-direction
  --Ly=<Ly>                               The spatial extent in the y-direction
  -h --help                               Show this screen.
"""


import pylab as pl
import pimchelp
from docopt import docopt
import loadgmt,kevent
import itertools

#cmap = loadgmt.getcmap('dca','alarm.p1.0.1')
#cmap = loadgmt.getcmap('cb/div','Spectral_11')
#cmap = loadgmt.getcmap('cb/div','RdGy_11')
cmap = loadgmt.getcmap('gist','earth')
cmap = loadgmt.getcmap('ncl','precip3_16lev')
#cmap = loadgmt.getcmap('ncl','temp_diff_18lev')

# ===================================================================
def main():

    # Read in the command line arguments
    args = docopt(__doc__)
    fileName = args['<input-file>']

    rskip = args['--rskip'] and float(args['--rskip'])

    # Open up the matrix file, and determine the size of the matrix
    data = pl.loadtxt(fileName,ndmin=2)
    N = int(pl.sqrt(data.shape[0]))

    # Get the spatial dimensions, otherwise just use N
    L = [N,N]
    for n,l in enumerate(['--Lx','--Ly']):
        if args[l]:
            L[n] = float(args[l])

    gridSize = []
    for cL in L:
        gridSize.append(cL/(1.0*N))

    # Reshape the data into a matrix
    density = data[:,1].reshape(N,N)
    ddensity = data[:,2].reshape(N,N)

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
            rho[n] += density[i,j]/(r*r)
            counts[n] += 1
        if rskip and r < rskip:
            density[i,j] = 0.0
    for n in range(num_rad_sep):
        if counts[n] > 0:
            rho[n] /= counts[n]

    # plot histograms in all three projections
    pl.figure(1)
    pl.imshow(density, cmap=cmap, 
              extent=[-0.5*L[0],0.5*L[0],-0.5*L[1],0.5*L[1]])
    pl.xlabel(r'$y\  [\AA]$')
    pl.ylabel(r'$x\  [\AA]$')
    pl.colorbar(shrink=0.4)

    r = pl.linspace(0,0.5*L[0],num_rad_sep)
    pl.figure(2)
    pl.plot(r,rho,'ro-', markersize=8)


    loc_ave_data = pl.zeros_like(density)
    for i,j in itertools.product(l,l):
        x = -0.5*L[0] + (i+0.5)*gridSize[0]
        y = -0.5*L[1] + (j+0.5)*gridSize[1]
        r = pl.sqrt(x*x + y*y)
        if r < 0.5*L[0]:
            count = 0
            for k,l in itertools.product(range(-1,2),range(-1,2)):
                if i+k < N and j+l < N:
                    loc_ave_data[i,j] += density[i+k,j+l]
                    count += 1
            if count > 0:
                loc_ave_data[i,j] /= count
        else:
            loc_ave_data[i,j] = density[i,j]

    pl.figure()
    pl.imshow(loc_ave_data, cmap=cmap, 
              extent=[-0.5*L[0],0.5*L[0],-0.5*L[1],0.5*L[1]])
    pl.xlabel(r'$y\  [\AA]$')
    pl.ylabel(r'$x\  [\AA]$')
    pl.colorbar(shrink=0.4)

    # Try to average neighboring bins together.
    ave = pl.zeros([N/2,N/2])
    for i,j in itertools.product(range(0,N-1,2),range(0,N-1,2)):
        for k,l in itertools.product(range(0,2),range(0,2)):
            ave[i/2,j/2] += density[i+k,j+l]
        ave[i/2,j/2] /= 4.0

    pl.figure()
    pl.imshow(ave, cmap=cmap, 
              extent=[-0.5*L[0],0.5*L[0],-0.5*L[1],0.5*L[1]])
    pl.xlabel(r'$y\  [\AA]$')
    pl.ylabel(r'$x\  [\AA]$')
    pl.colorbar(shrink=0.4)
   
    pl.show()
          
# ===================================================================
if __name__=="__main__":
    main()
