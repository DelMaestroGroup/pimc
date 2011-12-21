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
import plotoptions

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
    parser.add_argument('--plot','-p', help='The plot type. l = lines, p = points, \
                        e = errorbars', type=str, choices=['l','e','p','lp','le'], 
                        default='e')
    parser.add_argument('--ndim', '-d', help='Number of spatial dimensions.',
                        type=int, default=3) 
    parser.add_argument('--label', '-l', help='Parameter name for labels.', type=str)
    parser.add_argument('--xlim', '-x', help='x-axis limits', type=float,
                        nargs='+')
    parser.add_argument('--ylim', '-y', help='y-axis limits', type=float,
                        nargs='+')
    args = parser.parse_args()

    fileNames = args.fileNames
    estimatorToPlot = args.estimator
    NDIM = args.ndim
    plotType = args.plot
    varLabel = args.label
    xLim = args.xlim
    yLim = args.ylim

    if len(fileNames) < 1:
        parser.error("Need to specify at least one scalar estimator file")

    # Analyze the imput files
    reduce = pimchelp.ScalarReduce(fileNames,varLabel)

    # Get a color scheme and marker list
    numColors = max(reduce.getNumVarParams(),2)
    markers,colors  = plotoptions.markersColors(numColors)

    # create a description object
    descrip = pimchelp.Description(NDIM)

    # get the plot options
    pOptions = plotoptions.plotOptions(plotType)

    # Plot each estimator
    for n,est in enumerate(estimatorToPlot):
        pl.figure(n+1)
        ax = pl.subplot(111)
        for i in range(reduce.getNumVarParams()):
            lab = reduce.getVarLabel(i)
            pOptions['color'] = colors[i]
            if 'p' in plotType:
                pl.plot(reduce.param(), reduce.estimator(est,i),
                        markerfacecolor=colors[i],label=lab, **pOptions)
            if plotType == 'l':
                pl.plot(reduce.param(), reduce.estimator(est,i),
                        label=lab,**pOptions)
            elif 'e' in plotType:
                eb = pl.errorbar(reduce.param(), reduce.estimator(est,i),
                                 yerr=reduce.estimatorError(est,i), 
                                 markerfacecolor=colors[i], ecolor=colors[i],
                                 label=lab, **pOptions)
                
                # Set the width of the cap lines
                for cap in eb[1]:
                    cap.set_mew(2.0)

        #pl.tight_layout()
        pl.xlabel(descrip.paramLongName[reduce.reduceLabel])
        pl.ylabel(descrip.estimatorLongName[est])
        leg = pl.legend(frameon=False, loc='best', prop={'size':18})
        if xLim != None:
            pl.xlim(xLim[0],xLim[1])

        if yLim != None:
            pl.ylim(yLim[0],yLim[1])

    pl.show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

