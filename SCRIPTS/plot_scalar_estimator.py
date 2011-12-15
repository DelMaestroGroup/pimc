# plot_scalar_estimator.py
# Adrian Del Maestro
# 11.29.2011
# 
# Plot a single scalar estimator

import os,sys
import loadgmt,kevent
import argparse
import pyutils
import pimchelp
import numpy as np
import pylab as pl

# -----------------------------------------------------------------------------
def plotOptions(plotType):
    ''' Consctruct a dictionary of plot options based on what type of plot is
    desired.'''

    # First define the common properties
    pOps = {}

    # marker properties
    pOps['markersize'] = 9
    pOps['markeredgewidth'] = 1.0
    pOps['markeredgecolor'] = 'black'

    # line properties
    pOps['color'] = 'black'
    pOps['linewidth'] = 1.0

    if 'l' in plotType:
        pOps['linestyle'] = '-' 
    else:
        pOps['linestyle'] = 'None' 

    if plotType == 'p':
        pOps['linestyle'] = 'None'

    if plotType == 'l':
        pOps['linewidth'] = 2.0
        pOps['marker'] = None
        pOps['markeredgewidth'] = 0.0
        pOps['markersize'] = 0.0
        pOps.pop('color')

    if 'e' in plotType:
        pOps['capsize'] = 10
        pOps['elinewidth'] = 1.0
        
    return pOps


# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the argument parser
    parser = argparse.ArgumentParser(description='Plot a scalar estimator.')
    parser.add_argument('fileNames', help='Reduced scalar estimator files', 
                        nargs='+')
    parser.add_argument('--estimator','-e', help='A list of estimator names that \
                        are to be plotted.', type=str, nargs='+', required=True)
    parser.add_argument('--plot', help='The plot type. l = lines, p = points, \
                        e = errorbars', type=str, choices=['l','e','p','lp','le'], 
                        default='e')
    parser.add_argument('--ndim', '-d', help='Number of spatial dimensions.',
                        type=int, default=3) 
    args = parser.parse_args()

    fileNames = args.fileNames
    estimatorToPlot = args.estimator
    NDIM = args.ndim
    plotType = args.plot

    if len(fileNames) < 1:
        parser.error("Need to specify at least one scalar estimator file")

    # Analyze the imput files
    reduce = pimchelp.ScalarReduce(fileNames)

    # Get a color scheme and marker list
    # http://www.graphviz.org/content/color-names
    numColors = max(reduce.getNumVarParams(),5)
    colors  = loadgmt.getColorList('cb/div','Spectral_08',numColors)
#    colors.reverse()
    markers = loadgmt.getMarkerList()

    # create a description object
    descrip = pimchelp.Description(NDIM)

    # get the plot options
    pOptions = plotOptions(plotType)

    # Plot each estimator
    for n,est in enumerate(estimatorToPlot):
        pl.figure(n+1)
        ax = pl.subplot(111)
        for i in range(reduce.getNumVarParams()):
            lab = reduce.getVarLabel(i)
            if 'p' in plotType:
                pOptions['color'] = colors[i]
                pl.plot(reduce.param(), reduce.estimator(est,i),
                        marker=markers[0], markerfacecolor=colors[i],label=lab,
                        **pOptions)
            if plotType == 'l':
                pl.plot(reduce.param(), reduce.estimator(est,i),
                        color=colors[i], label=lab,**pOptions)
            elif 'e' in plotType:
                eb = pl.errorbar(reduce.param(), reduce.estimator(est,i),
                                 yerr=reduce.estimatorError(est,i), 
                                 marker=markers[i], markerfacecolor=colors[i],
                                 ecolor=colors[i], label=lab, **pOptions)
                
                # Set the width of the cap lines
                for cap in eb[1]:
                    cap.set_mew(1.0)

        #pl.tight_layout()
        pl.xlabel(descrip.paramLongName[reduce.reduceLabel])
        pl.ylabel(descrip.estimatorLongName[est])
        leg = pl.legend(frameon=False, loc='best', prop={'size':18})

    pl.show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

