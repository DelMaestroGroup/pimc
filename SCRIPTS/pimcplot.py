#!/usr/bin/python
# pimcplot.py
# Adrian Del Maestro
# 07.20.2009
# 
# Plot rough estimators vs. MC steps for files supplied as input

import matplotlib
matplotlib.use('TKAgg')

import os,sys
import pyutils
from optparse import OptionParser
import loadgmt,kevent
from pylab import *


# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = OptionParser() 
    parser.add_option("-c", "--column", dest="col", type="int",\
            help="which column should we plot?")

    # parse the command line options and get the file name
    (options, args) = parser.parse_args() 
    if len(args) < 1: 
        parser.error("need a file name")

    if (not options.col):
        parser.error("need to supply a column to plot!")

    options.col -= 1
    col = [options.col]
    
    fileNames = args
    numFiles = len(fileNames)

    # We count the number of lines in the estimator file to make sure we have
    # some data and grab the headers
    estFile = open(fileNames[0],'r');
    estLines = estFile.readlines();
    numLines = len(estLines) - 2    # We expect two comment lines
    pimcid = estLines[0]
    headers = estLines[1].split()
    estFile.close()

    # ============================================================================
    # Figure 1 : column vs. MC Steps
    # ============================================================================
    figure(1)
    connect('key_press_event',kevent.press)

    colors  = loadgmt.getColorList('cw/1','cw1-029',max(numFiles,2))

    n = 0
    for fileName in fileNames:

        dataFile = open(fileName,'r');
        dataLines = dataFile.readlines();

        if len(dataLines) > 2:
            data = loadtxt(fileName,usecols=col)
            plot(data,marker='s',color=colors[n],markeredgecolor=colors[n],\
                        markersize=4,linestyle='-',linewidth=1.0)
    
            n += 1

    ylabel(r'$%s$' % headers[options.col+1])
    xlabel("MC Steps")

    # ============================================================================
    # Figure 2 : running average of column vs. MC Steps
    # ============================================================================
    figure(2)
    connect('key_press_event',kevent.press)

    n = 0
    for fileName in fileNames:

        dataFile = open(fileName,'r');
        dataLines = dataFile.readlines();

        if len(dataLines) > 2:

            data = loadtxt(fileName,usecols=col)
            if size(data) > 1:
                aveData = []
                for r in range(1,len(data)+1):
                    aveData.append(mean(data[0:r]))

                plot(aveData,color=colors[n],linewidth=1.5,marker='None',linestyle='-')
    
                n += 1

    ylabel(r'$\langle %s \rangle$' % headers[options.col+1])
    xlabel("MC Steps")
    tight_layout()

    show()
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
