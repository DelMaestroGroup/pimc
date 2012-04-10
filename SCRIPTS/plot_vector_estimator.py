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
import plotoptions

# --------------------------------------------------------------------------------
def getIndex(A,N,NDIM):
    ''' Returns the reconstructed index of a 2d array'''
    index = []
    for i in range(NDIM):
        fact = 1
        for j in range(i+1,NDIM):
            fact *= N[j]
        index.append((A/fact) % N[i])
    return index

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the argument parser
    parser = argparse.ArgumentParser(description='Plot a vector estimator.')
    parser.add_argument('fileNames', help='Reduced scalar estimator files', nargs='+')
    parser.add_argument('--ndim', '-d', help='Number of spatial dimensions.',
                        type=int, default=3) 
    parser.add_argument('--plot','-p', help='The plot type. l = lines, p = points, \
                        e = errorbars', type=str, choices=['l','e','p','lp','le'], 
                        default='e')
    parser.add_argument('--label', '-l', help='Parameter name for labels.', type=str)
    parser.add_argument('--subplot', '-s', help='Dimensions of a subplot matrix.', 
                        type=int, nargs='+')
    parser.add_argument('--xlim', '-x', help='x-axis limits', type=float,
                        nargs='+')
    parser.add_argument('--ylim', '-y', help='y-axis limits', type=float,
                        nargs='+')
    parser.add_argument('--log', help='Set y-axis log scale',
                        action='store_true')
    args = parser.parse_args()

    fileNames = args.fileNames
    NDIM = args.ndim
    varLabel = args.label
    plotType = args.plot
    subplot = args.subplot
    xLim = args.xlim
    yLim = args.ylim
    log = args.log

    if len(fileNames) < 1:
        parser.error("Need to specify at least one vector estimator file")

    # Analyze the imput files
    estimatorName = pimchelp.getVectorEstimatorName(fileNames[0])
    reduce = pimchelp.VectorReduce(fileNames,estimatorName,varLabel)

    # Get a color scheme and marker list
    numColors = reduce.getNumReduceParams()
    markers,colors  = plotoptions.markersColors(numColors)

    # get the plot options
    pOptions = plotoptions.plotOptions(plotType)

    # create a description object
    descrip = pimchelp.Description(NDIM)

    # Plot each estimator
    for varIndex in range(reduce.getNumVarParams()):
        pl.figure(varIndex+1)
        #ax = pl.subplot(111)

        for reduceIndex in range(reduce.getNumReduceParams()):
            lab = reduce.getReduceLabel(reduceIndex)
            pOptions['color'] = colors[reduceIndex]
            if 'e' in plotType:
                eb = pl.errorbar(reduce.x(varIndex,reduceIndex),
                                 reduce.estimator(varIndex,reduceIndex),
                                 yerr=reduce.estimatorError(varIndex,reduceIndex), 
                                 markerfacecolor=colors[reduceIndex],
                                 ecolor=colors[reduceIndex], label=lab, 
                                 **pOptions)
                
                # Set the width of the cap lines
                for cap in eb[1]:
                    cap.set_mew(1.0)
            else:
                pl.plot(reduce.x(varIndex,reduceIndex),
                        reduce.estimator(varIndex,reduceIndex),
                        label=lab, **pOptions)


        pl.xlabel(descrip.estimatorXLongName[estimatorName])
        pl.ylabel(descrip.estimatorLongName[estimatorName])
        leg = pl.legend(frameon=False, loc='best', prop={'size':18})

#        # Add a colorbar
#        cmap = loadgmt.getcmap('cb/div','Spectral_08')
#        T = reduce.param()[::2]
#        cb = pl.colorbar(loadgmt.getMap(cmap,T),ax=ax,ticks=T)
#        cb.ax.set_ylabel('Temperature',rotation=-90)

        if xLim != None:
            pl.xlim(xLim[0],xLim[1])

        if yLim != None:
            pl.ylim(yLim[0],yLim[1])

        # Set a log scale
        if log:
            ax = pl.subplot(111)
            ax.set_yscale('log')


    # Plot a possible subplot matrix
    if subplot != None:
        f, ax = pl.subplots(subplot[0], subplot[1], sharex=True, squeeze=False, sharey=True)
        numReduce = reduce.getNumReduceParams()

        for reduceIndex in range(reduce.getNumReduceParams()):
            id = getIndex(reduceIndex,subplot,2)
            pOptions['color'] = colors[reduceIndex]
            lab = reduce.getReduceLabel(reduceIndex)

            for varIndex in range(reduce.getNumVarParams()):
                pOptions['marker'] = markers[varIndex]

                if 'e' in plotType:
                    eb = ax[id[0],id[1]].errorbar(reduce.x(varIndex,reduceIndex),
                                                  reduce.estimator(varIndex,reduceIndex),
                                                  yerr=reduce.estimatorError(0,reduceIndex), 
                                                  markerfacecolor=colors[reduceIndex],
                                                  ecolor=colors[reduceIndex], label=lab, 
                                                  **pOptions)
                    
                    # Set the width of the cap lines
                    for cap in eb[1]:
                        cap.set_mew(1.0)
                else:
                    ax[id[0],id[1]].plot(reduce.x(varIndex,reduceIndex),
                                         reduce.estimator(varIndex,reduceIndex),
                                         label=lab, **pOptions)

            ax[id[0],id[1]].legend(frameon=False, loc='upper left', prop={'size':18})
            if xLim != None:
                ax[id[0],id[1]].set_xlim(xLim[0],xLim[1])

            if yLim != None:
                ax[id[0],id[1]].set_ylim(yLim[0],yLim[1])

        [ax[-1,n].set_xlabel(descrip.estimatorXLongName[estimatorName]) for n in range(subplot[1])]
        [ax[n,0].set_ylabel(descrip.estimatorLongName[estimatorName]) for n in range(subplot[0])]
        f.subplots_adjust(hspace=0.10)
        f.subplots_adjust(wspace=0.10)

    pl.show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

