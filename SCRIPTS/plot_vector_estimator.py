# plot_vector_estimator.py
# Adrian Del Maestro
# 12.02.2011
# 
# Plot possibly many vector estimators

import os,sys
import loadgmt,kevent
import argparse
import pyutils
import pimchelp
import numpy as np
import pylab as pl

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the argument parser
    parser = argparse.ArgumentParser(description='Plot a vector estimator.')
    parser.add_argument('fileNames', help='Reduced scalar estimator files', nargs='+')
    parser.add_argument('--estimator','-e', help='A list of estimator names that \
                        are to be plotted.', type=str, nargs='+')
    parser.add_argument('--no_error', help='Turn off error bars.', action='store_true')
    parser.add_argument('--lines', help='Plot with lines.', action='store_true')
    parser.add_argument('--ndim', '-d', help='Number of spatial dimensions.',
                        type=int, default=3) 
    args = parser.parse_args()

    fileNames = args.fileNames
    estimatorToPlot = args.estimator
    noErrorBars = args.no_error
    NDIM = args.ndim
    lines = args.lines

    if len(fileNames) < 1:
        parser.error("Need to specify at least one vector estimator file")

    # Analyze the imput files
    estimatorName = pimchelp.getVectorEstimatorName(fileNames[0])
    reduce = pimchelp.VectorReduce(fileNames,estimatorName)

    # Get a color scheme and marker list
    # http://www.graphviz.org/content/color-names
    numColors = max(reduce.getNumReduceParams(),5)
    colors  = loadgmt.getColorList('cb/div','Spectral_05',numColors)
#    colors.reverse()
    markers = loadgmt.getMarkerList()

    # create a description object
    descrip = pimchelp.Description(NDIM)

    # Plot each estimator
    for varIndex in range(reduce.getNumVarParams()):
        pl.figure(varIndex+1)
        for reduceIndex in range(reduce.getNumReduceParams()):
            lab = reduce.getReduceLabel(reduceIndex)
            if lines:
                pl.plot(reduce.x(varIndex,reduceIndex),
                        reduce.estimator(varIndex,reduceIndex),
                        marker=None, label=lab, linestyle='-',
                        markeredgewidth=0,markersize=0,
                        color=colors[reduceIndex], linewidth=2.0)

            else:
                eb = pl.errorbar(reduce.x(varIndex,reduceIndex),
                                 reduce.estimator(varIndex,reduceIndex),
                                 yerr=reduce.estimatorError(varIndex,reduceIndex), 
                                 marker=markers[0], markersize=10,
                                 markerfacecolor=colors[reduceIndex],
                                 markeredgewidth=0.2, markeredgecolor='k', capsize=10, 
                                 ecolor=colors[reduceIndex], elinewidth=1.0, label=lab, 
                                 linestyle='-', color='k', linewidth=0.2)
                
                # Set the width of the cap lines
                for cap in eb[1]:
                    cap.set_mew(1.0)


        #pl.tight_layout()
        pl.xlabel(descrip.estimatorXLongName[estimatorName])
        pl.ylabel(descrip.estimatorLongName[estimatorName])
        pl.legend(frameon=False, loc='best')
    pl.show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

