# reduce-one.py
# Adrian Del Maestro
# 09.03.2009
# 
# Reduce and average results for a single PIMC run based on a single parameter
# which varies.  This could be density for fixed system size, etc.

import os,sys,glob
import loadgmt,kevent
import pimchelp
from optparse import OptionParser
from pylab import *

# ----------------------------------------------------------------------
def getStats(data,dim=0):
    ''' Get the average and error of all columns in the data matrix. '''

    if ndim(data) > dim:
        numBins  = size(data,dim) 
        dataAve  = average(data,dim) 
        dataAve2 = average(data*data,dim) 
        dataErr   = sqrt( abs(dataAve2-dataAve**2)/(1.0*numBins-1.0) ) 
    else:
        dataAve = data
        dataErr = 0.0*data

    return dataAve,dataErr

# -----------------------------------------------------------------------------
def getScalarEst(type,pimc,outName,reduceFlag):
    ''' Return the arrays containing the reduced averaged scalar
        estimators in question.'''

    fileNames = pimc.getFileList(type)
    headers   = pimchelp.getHeadersFromFile(fileNames[0])

    ave = zeros([len(fileNames),len(headers)],float)
    err = zeros([len(fileNames),len(headers)],float)
    for i,fname in enumerate(fileNames):
        # Compute the averages and error
        data = loadtxt(fname)
        ave[i,:],err[i,:] = getStats(data)
    
    # output the estimator data to disk
    outFile = open('%s-%s' % (type,outName),'w');

    # the headers
    outFile.write('#%15s' % reduceFlag[0])
    for head in headers:
        outFile.write('%16s%16s' % (head,'+/-'))
    outFile.write('\n')

    # the data
    for i,f in enumerate(fileNames):
        outFile.write('%16.8E' % float(pimc.params[pimc.id[i]][reduceFlag[1]]))
        for j,h in enumerate(headers):
            outFile.write('%16.8E%16.8E' % (ave[i,j],err[i,j]))
        outFile.write('\n')
    outFile.close()

    return headers,ave,err;

# -----------------------------------------------------------------------------
def getVectorEst(type,pimc,outName,reduceFlag,xlab,ylab):
    ''' Return the arrays consisting of the reduec averaged vector 
        estimators. '''

    fileNames = pimc.getFileList(type)
    headers   = pimchelp.getHeadersFromFile(fileNames[0])

    numParams = len(fileNames)
    Nx = len(headers)

    x   = zeros([numParams,Nx],float)
    ave = zeros([numParams,Nx],float)
    err = zeros([numParams,Nx],float)

    
    for i,fname in enumerate(fileNames):

        # Get the estimator data and compute averages
        data = loadtxt(fname)
        ave[i,:],err[i,:] = getStats(data)

        # get the headers
        x[i,:] = pimchelp.getHeadersFromFile(fname)

        # Compute the normalized averages and error for the OBDM
        if type == 'obdm':
            norm = ave[i,0]
            ave[i,:] /= norm
            err[i,:] /= norm

    # output the vector data to disk
    outFile = open('%s-%s' % (type,outName),'w');

    # the headers
    lab = '%s = %4.2f' % (reduceFlag[0],float(pimc.params[pimc.id[0]][reduceFlag[1]]))
    outFile.write('#%15s%16s%16s' % ('',lab,''))
    for j in range(numParams-1):
        lab = '%s = %4.2f' % (reduceFlag[0],float(pimc.params[pimc.id[j+1]][reduceFlag[1]]))
        outFile.write('%16s%16s%16s' % ('',lab,''))
    outFile.write('\n')
    outFile.write('#%15s%16s%16s' % (xlab,ylab,'+/-'))
    for j in range(numParams-1):
        outFile.write('%16s%16s%16s' % (xlab,ylab,'+/-'))
    outFile.write('\n')

    # the data
    for i,h in enumerate(headers):
        for j in range(numParams):
            outFile.write('%16.8E%16.8E%16.8E' % (x[j,i],ave[j,i],err[j,i]))
        outFile.write('\n')
    outFile.close()

    return x,ave,err

# -----------------------------------------------------------------------------
def getKappa(pimc,outName,reduceFlag):
    ''' Return the arrays containing the reduced averaged compressibility. '''

    fileNames = pimc.getFileList('estimator')
    headers   = pimchelp.getHeadersFromFile(fileNames[0])

    aveKappa = zeros([len(fileNames)],float)
    errKappa = zeros([len(fileNames)],float)

    for i,fname in enumerate(fileNames):

        # Now get the temperature, volume and linear system size
        ID = pimc.getID(fname)
        T = float(pimc.params[ID]['Temperature'])
        V = float(pimc.params[ID]['Container Volume'])
        box = pimc.params[ID]['Container Dimensions'].split('x')
        L = float(box[-1])

        # Compute the average compressibility and its error
        estData = loadtxt(fname)

        N     = estData[:,headers.index('N')]
        N2    = estData[:,headers.index('N^2')] 

        numBins = len(N)
            
        # Get the averages
        aveN,errN = getStats(N)
        aveN2,errN2 = getStats(N2)
        aveNN2,errNN2 = getStats(N*N2)

        # Get the covariance
        # This is finite, so it must be calculated!
        covNN2 = (aveNN2 - aveN*aveN2)/(1.0*numBins-1.0)

        # Get the value of rho^2 * kappa and the error
        aveKappa[i] = (aveN2-aveN**2)/(T*V)
        errKappa[i] = sqrt(errN2**2 + 4.0*errN**2*aveN**2 - 4.0*aveN*covNN2)/(T*V)
    
    # output the estimator data to disk
    outFile = open('%s-%s' % ('kappa',outName),'w');

    # the headers
    outFile.write('#%15s' % reduceFlag[0])
    outFile.write('%16s%16s' % ('kappa','+/-'))
    outFile.write('\n')

    # the data
    for i in range(len(fileNames)):
        outFile.write('%16.8E' % float(pimc.params[pimc.id[i]][reduceFlag[1]]))
        outFile.write('%16.8E%16.8E\n' % (aveKappa[i],errKappa[i]))
    outFile.close()

    return aveKappa,errKappa


# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # define the mapping between short names and label names 
    shortFlags = ['n','T','N','t','u','V']
    parMap = {'n':'Initial Density', 'T':'Temperature', 'N':'Initial Number Particles',\
            't':'Imaginary Time Step', 'u':'Chemical Potential', 'V':'Container Volume'}

    # setup the command line parser options 
    parser = OptionParser() 
    parser.add_option("-T", "--temperature", dest="T", type="float", \
            help="simulation temperature in Kelvin") 
    parser.add_option("-N", "--number-particles", dest="N", type="int",\
            help="number of particles") 
    parser.add_option("-n", "--density", dest="n", type="float",\
            help="number density in Angstroms^{-d}")
    parser.add_option("-t", "--imag-time-step", dest="tau", type="float",\
            help="imaginary time step")
    parser.add_option("-u", "--chemical-potential", dest="mu", type="float",\
            help="chemical potential in Kelvin") 
    parser.add_option("-V", "--volume", dest="V", type="float",\
            help="volume in Angstroms^d") 
    parser.add_option("-r", "--reduce", dest="reduce", choices=['T','N','n','u','t','V'],\
            help="variable name for reduction [T,N,n,u,M,V]") 
    parser.add_option("-g", "--grand-canonical", action="store_true", dest="gce",\
            help="are we in the grand canonical ensemble?")
    parser.add_option("-o", "--obdm", action="store_false", dest="obdm",\
            help="do we want to measure the obdm?") 
    parser.add_option("-p", "--plot", action="store_true", dest="plot",\
            help="do we want to produce data plots?") 
    parser.set_defaults(gce=False)
    parser.set_defaults(plot=False)
    parser.set_defaults(obdm=True)

    # parse the command line options and get the reduce flag
    (options, args) = parser.parse_args() 
    if len(args) > 0: 
        parser.error("incorrect number of arguments")
    
    if (not options.reduce):
        parser.error("need a correct reduce flag (-r,--reduce): [T,N,n,u,t,V]")

    # Check that we are in the correct ensemble
    pimchelp.checkEnsemble(options.gce)

    dataName,outName = pimchelp.getFileString(options)
    reduceFlag = []
    reduceFlag.append(options.reduce)
    reduceFlag.append(parMap[options.reduce])

    # Create the PIMC analysis helper and fill up the simulation parameters maps
    pimc = pimchelp.PimcHelp(dataName,options.gce)
    pimc.getSimulationParameters()

    # Form the full output file name
    outName += '.dat'

    # We first reduce the scalar estimators and output them to disk
    head1,scAve1,scErr1 = getScalarEst('estimator',pimc,outName,reduceFlag)
    head2,scAve2,scErr2 = getScalarEst('super',pimc,outName,reduceFlag)

    # Now we do the normalized one body density matrix
    if options.obdm:
        x1,ave1,err1 = getVectorEst('obdm',pimc,outName,reduceFlag,'r [A]','n(r)')

    # Now we do the pair correlation function
    x2,ave2,err2 = getVectorEst('pair',pimc,outName,reduceFlag,'r [A]','g(r)')

    # If we are reducing for the case of a cylindrical geometry
    if len(glob.glob('CYLINDER')) > 0:
        x3,ave3,err3 = getVectorEst('radial',pimc,outName,reduceFlag,'r [A]','rho(r)')

    # Compute the number distribution function and compressibility if we are in
    # the grand canonical ensemble
    if options.gce:
        x4,ave4,err4 = getVectorEst('number',pimc,outName,reduceFlag,'N','P(N)')
        kappa,kappaErr = getKappa(pimc,outName,reduceFlag)

    # Do we show plots?
    if options.plot:

        # -----------------------------------------------------------------------------
        # Plot the averaged data
        # -----------------------------------------------------------------------------
    
        # Get the changing parameter that we are plotting against
        param = []
        for ID in pimc.id:
            param.append(float(pimc.params[ID][reduceFlag[1]]))
    
        numParams = len(param)
        markers = loadgmt.getMarkerList()
        headLab = ['E/N','K/N','V/N','N', 'diagonal']
        dataCol = []
        for head in headLab:
            n = 0
            for h in head1:
                if head == h:
                    dataCol.append(n)
                    break
                n += 1
        yLabelCol = ['Energy / N', 'Kinetic Energy / N', 'Potential Energy / N',\
                'Number Particles', 'Diagonal Fraction']

        colors  = loadgmt.getColorList('cw/1','cw1-029',max(numParams,len(headLab)))
    
        # ============================================================================
        # Figure -- Various thermodynamic quantities
        # ============================================================================
        figNum = 1
        for n in range(len(dataCol)):
            figure(figNum)
            connect('key_press_event',kevent.press)
    
            errorbar(param, scAve1[:,dataCol[n]], yerr=scErr1[:,dataCol[n]],\
                    color=colors[n],marker=markers[n],markeredgecolor=colors[n],\
                    markersize=8,linestyle='None',capsize=4)
    
            xlabel('%s'%options.reduce)
            ylabel(yLabelCol[n])
            tight_layout()
            #savefig('tba-energy.eps')
            figNum += 1
    
        # ============================================================================
        # Figure -- The superfluid density
        # ============================================================================
        figure(figNum)
        connect('key_press_event',kevent.press)
    
        errorbar(param, scAve2[:,0], yerr=scErr2[:,0],\
                color=colors[0],marker=markers[0],markeredgecolor=colors[0],\
                markersize=8,linestyle='None',capsize=4)
    
        tight_layout()
        xlabel('%s'%options.reduce)
        ylabel('Superfluid Density')
    
        # ============================================================================
        # Figure -- The one body density matrix
        # ============================================================================
        if options.obdm:
            figNum += 1
            figure(figNum)
            connect('key_press_event',kevent.press)
            ax = subplot(111)
    
            for n in range(numParams):
                lab = '%s = %s' % (options.reduce,param[n])
                errorbar(x1[n,:], (ave1[n,:]+1.0E-15), err1[n,:],color=colors[n],marker=markers[n],\
                        markeredgecolor=colors[n], markersize=8,linestyle='None',label=lab)
    
                #axis([0,21,1.0E-5,1.1])
            xlabel('r [Angstroms]')
            ylabel('One Body Density Matrix')
            tight_layout()
            legend(loc='best', frameon=False, prop={'size':16},ncol=2)
            #savefig('tba-obdm.eps')
    
        # ============================================================================
        # Figure -- The pair correlation function
        # ============================================================================
        figNum += 1
        figure(figNum)
        connect('key_press_event',kevent.press)
    
        for n in range(numParams):
            lab = '%s = %s' % (options.reduce,param[n])
            errorbar(x2[n,:], ave2[n,:], yerr=err2[n,:],color=colors[n],marker=markers[0],\
                    markeredgecolor=colors[n], markersize=8,linestyle='None',label=lab,capsize=6)
    
            #   axis([0,256,1.0E-5,1.2])
        xlabel('r [Angstroms]')
        ylabel('Pair Correlation Function')
        legend(loc='best', frameon=False, prop={'size':16},ncol=2)
        tight_layout()
        #savefig('tba-pair.eps')
    
        # We only plot the compressibility if we are in the grand-canonical ensemble
        if options.gce: 
    
            # ============================================================================
            # Figure -- The Number distribution
            # ============================================================================
            figNum += 1
            figure(figNum)
            connect('key_press_event',kevent.press) 

            # Find which column contains the average number of particles
            for hn,h in enumerate(head1):
                if h == 'N':
                    break

            for n in range(numParams): 
                lab = '%s = %s' % (options.reduce,param[n]) 
                aN = scAve1[n,hn] 
                errorbar(x4[n,:]-aN, ave4[n,:], err4[n,:],color=colors[n],marker=markers[0],\
                         markeredgecolor=colors[n],\
                         markersize=8,linestyle='None',label=lab,capsize=6) 
    
            axis([-30,30,0.0,1.2])
            xlabel('N-<N>')
            ylabel('P(N)')
            tight_layout()
            legend(loc='best', frameon=False, prop={'size':16},ncol=2)
    
            # ============================================================================
            # Figure -- The Compressibility
            # ============================================================================
            figNum += 1
            figure(figNum)
            connect('key_press_event',kevent.press)

            errorbar(param, kappa, yerr=kappaErr, color=colors[0],marker=markers[0],\
                    markeredgecolor=colors[0], markersize=8,linestyle='None',capsize=6)
    
            tight_layout()
            xlabel('%s'%options.reduce)
            ylabel(r'$\rho^2 \kappa$')
    
        # ============================================================================
        # Figure -- The radial density
        # ============================================================================
        if len(glob.glob('CYLINDER')) > 0:
            figNum += 1
            figure(figNum)
            connect('key_press_event',kevent.press)
            ax = subplot(111)
    
            for n in range(numParams):
                lab = '%s = %s' % (options.reduce,param[n])
                errorbar(x3[n,:], (ave3[n,:]+1.0E-15), err3[n,:],color=colors[n],marker=markers[0],\
                        markeredgecolor=colors[n], markersize=8,linestyle='None',label=lab)
    
                #axis([0,21,1.0E-5,1.1])
            tight_layout()
            xlabel('r [Angstroms]')
            ylabel('Radial Density')
            legend(loc='best', frameon=False, prop={'size':16},ncol=2)
    
        show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

