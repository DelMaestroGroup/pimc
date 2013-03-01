''' PimcHelp - helper methods for analzing pimc output data.
'''
import os
import glob
from operator import itemgetter, attrgetter

# ----------------------------------------------------------------------
def getVectorEstimatorName(fileName):
    '''Determine the name of a reduced estimator.'''

    # Get the actual file name and then split at the '-' returning the first
    # word
    if '/' in fileName:
        estName = fileName.split('/')[-1].split('-')[0]
    else:
        estName = fileName.split('-')[0]

    return estName

# ----------------------------------------------------------------------------
def getFileNameParameters(fname):
    '''Get the parameters from the output filename.  
    
    We need to be careful due to the problems associated with splitting at the
    '-' character when there are potential minus signs.
    '''

    fileParts  = fname.rstrip('.dat').split('-')
    pIndex = []
    for n,part in enumerate(fileParts):
        if part == '':
            fileParts[n+1] = '-' + fileParts[n+1]
            pIndex.append(n)
    for n in pIndex:
        fileParts.pop(n)

    return fileParts

# -----------------------------------------------------------------------------
# !!! BROKEN FOR CANONICAL ENSEMBLE
# !!! due to different file labelling scheme
# -----------------------------------------------------------------------------
def sortFileNames(fileNames): 
    '''Try to sort filenames by the numeric values of their parameter
    strings.'''

    fileTuples = []
    for fname in fileNames:
        # break up the parameters in the file name
        fileParts = getFileNameParameters(fname)

        # get the tuple
        tup = (fname, fileParts[0], fileParts[1], float(fileParts[2]),
               float(fileParts[3]), float(fileParts[4]), float(fileParts[5]),
               int(fileParts[6]))
        fileTuples.append(tup)

    # sort by keys for each column from left to right
    for n in range(1,7):
        fileTuples = sorted(fileTuples, key=itemgetter(n))
    
    # this is the new fangled way, but it is not backwards compatible
    #fileTuples = sorted(fileTuples, key=itemgetter(1,2,3,4,5,6))
    
    # get the sorted file names
    sortedFileNames = []
    for ft in fileTuples:
        sortedFileNames.append(ft[0])


    # return the sorted file names
    return sortedFileNames

# -----------------------------------------------------------------------------
def getParFromReduceFile(fileName): 
    ''' Create a parameter map out of an output file name 
        (coming from reduce-one.py)'''

    # We split the fileName at reduce, and grab the latter half
    dataName = fileName.partition('reduce-')[2].rstrip('.dat').split('-')

    # This part fixes any difficulties related to a negative chemical
    # potential
    for data in dataName:
        if data == '':
            n = dataName.index(data)
            dataName.pop(n)
            dataName[n] = '-' + dataName[n]

    # Now, parse the file name and get the parameter values
    dataMap = {}
    n = 0
    while n < len(dataName):
        if dataName[n] == 'N':
            dataMap[dataName[n]]= int(dataName[n+1])
        else:
            dataMap[dataName[n]]= float(dataName[n+1])
        n += 2

    return dataMap

# -----------------------------------------------------------------------------
def getParFromPIMCFile(fileName): 
    ''' Create a parameter map out of an output file name of the form
    gce-*-T-V-mu-tau-id.dat or ce-*-T-N-n-tau-id.dat.'''

    # We split the fileName at reduce, and grab the latter half
    dataName = fileName.strip('.dat').split('-')
    for data in dataName:
        if data == '':
            n = dataName.index(data)
            dataName.pop(n)
            dataName[n] = '-' + dataName[n]

    # Now, parse the file name and get the parameter values
    dataMap = {}

    dataMap['T'] = float(dataName[2])
    dataMap['tau'] = float(dataName[5])
    dataMap['id'] = int(dataName[6])

    # The volume is either in the file name, or computed depending
    # on whether we are in the GCE
    if 'gce' in dataName[0]:
        dataMap['V'] = float(dataName[3])
        dataMap['mu'] = float(dataName[4])
        dataMap['gce'] = True
    else:
        dataMap['V'] = float(dataName[3])/float(dataName[4])
        dataMap['N'] = int(dataName[3])
        dataMap['n'] = float(dataName[4])
        dataMap['gce'] = False

    return dataMap

# -----------------------------------------------------------------------------
def getHeadersFromFile(fileName, skipLines=0): 
    ''' Get the data column headers from a PIMC output file. '''

    with open(fileName,'r') as inFile:
        inLines = inFile.readlines();
        n = skipLines
        if 'PIMCID' in inLines[n]:
            headers = inLines[n+1].split()
        else:
            headers = inLines[n].split()
        headers.pop(0)

    return headers

# -----------------------------------------------------------------------------
def getHeadersDict(fileName, removeLab=None, skipLines=0): 
    '''Construct a header dictionary from a filename.'''

    # Get all headers
    headers = getHeadersFromFile(fileName,skipLines)

    if removeLab != None:
        headers.remove(removeLab)

    headDict = {}
    for n,head in enumerate(headers):
        if head != '+/-':
            headDict[head] = n
            if n < (len(headers) - 1) and headers[n+1] == '+/-':
                headDict['d_' + head] = n+1

    return headDict

# -------------------------------------------------------------------------------
def checkEnsemble(canonical):
    ''' Here we make sure that the correct ensemble flag is specified. '''

    import sys
    gceFiles = (len(glob.glob('gce-*'))>0)
    ceFiles = (len(glob.glob('ce-*'))>0)

    if (ceFiles and not gceFiles) and not canonical:
        sys.exit('Need to include --canonical for the canonical ensemble!')

# -------------------------------------------------------------------------------
def getFileString(options,reduce=True):
    ''' Using the command line flags, form the input file string that will
        be used to open all data files. '''

    # we simply go through all possible options and figure out what the
    # reduce variable is
    out = ""
    if options.T is not None:
        flagT = "%06.3f" % options.T
        out += '-T-%s' % flagT
    else:
        flagT = "*"

    if options.N is not None:
        flagN = "%04d" % options.N
        out += '-N-%s' % flagN
    else:
        flagN = "*"

    if options.n is not None:
        flagn = "%06.3f" % options.n
        out += '-n-%s' % flagn
    else:
        flagn = "*"

    if options.tau is not None:
        flagtau = "%7.5f" % options.tau
        out += '-t-%s' % flagtau
    else:
        flagtau = "*"

    if options.mu is not None:
        flagmu = "%+08.3f" % options.mu
        out += '-u-%s' % flagmu
    else:
        flagmu = "*"

    if options.L is not None:
        flagL = "%07.3f" % options.L
        out += '-L-%s' % flagL
    else:
        flagL = "*"

    if options.canonical:
        dataName = '%s-%s-%s-%s-*.dat' % (flagT,flagN,flagn,flagtau)
    else:
        dataName = '%s-%s-%s-%s-*.dat' % (flagT,flagL,flagmu,flagtau)

    if reduce:
        outName = '%s-reduce%s' % (options.reduce,out)
        return dataName,outName
    else:
        return dataName

# -------------------------------------------------------------------------------
class PimcHelp:
    ''' Helper methods for analzing pimc output data. '''

    # ----------------------------------------------------------------------
    def __init__(self,dataName,canonical,baseDir=''):

        self.dataName = dataName 
        self.baseDir  = baseDir
        self.params = {}
        self.id = []
        if canonical:
            self.prefix='ce'
        else:
            self.prefix='gce'

        # The data file and all output file names
        self.dataType = ['estimator', 'obdm', 'pair', 'pcycle', 'super', 'worm', 
                         'radial', 'radwind', 'radarea', 'planedensity',
                         'planewind', 'planearea']
        if not canonical:
            self.dataType.append('number')
        self.outType  = ['estimator', 'number', 'obdm', 'pair', 'pcycle', 'super', 
                         'worm', 'radial', 'log', 'state']

    # -----------------------------------------------------------------------------
    def getID(self,fileName): 
        ''' Return the ID number corresponding to a given filename. '''
        #ID = int(fileName.rstrip('.dat').split('-')[-1])
        ID = int(fileName[-13:-4])
        return ID

    # ----------------------------------------------------------------------
    def getParameterMap(self,logName): 
        '''Given a log file name, return the parameter map. '''

        # Get the values of all simulation parameters
        paramsMap = {}
        params = False
        with open(logName, 'r') as logFile:
            for line in logFile:
                if 'Begin Simulation Parameters' in line:
                    params = True
                elif 'End Simulation Parameters' in line:
                    break

                if params and ':' in line:
                    keyVal = line.split(':')
                    paramsMap[keyVal[0].strip()] = keyVal[1].strip()

        # Add an element to the parameter map for the linear dimension (Lz) of
        # the container
        paramsMap['Container Length'] = paramsMap['Container Dimensions'].split('x')[-1]

        return paramsMap

    # -----------------------------------------------------------------------------
    def getSimulationParameters(self): 
        '''Get the full list of parameter maps for each input file and a list of
           ID numbers. '''

        # Get the list of log files
        fileNames = self.getFileList("log")

        self.params = {}
        self.id = []
        for fname in fileNames:
            ID = self.getID(fname)
            self.id.append(ID)
            self.params[ID] = self.getParameterMap(fname)

    # -----------------------------------------------------------------------------
    # def getFileInfo(self,type):
    #     ''' Get the names of the input files, how many of them there are, 
    #         and how many columns of data they contain.'''

    #     fileNames = self.getFileList(type)

    #     # The number of parameter files
    #     numParams = len(fileNames);

    #     # Open a sample data file and count the number of columns
    #     inFile = open(fileNames[0],'r')
    #     lines = inFile.readlines();
    #     for line in lines:
    #         if not (line[0] == '#'):
    #             cols = line.split()
    #             break
    #     numDataCols = len(cols)
    #     
    #     return fileNames,numParams,numDataCols

    # -----------------------------------------------------------------------------
    def getFileList(self,type,idList=None):
        ''' Get a list of input files based on their type, and possibly a number
            of unique ID's'''

        fileNames = []

        # We want all the file names here
        if not idList:
            fileLoc = '%s%s-%s-%s' % (self.baseDir,self.prefix,type,self.dataName)
            fileNames = glob.glob(fileLoc)

            # Now sort them
            fileNames  = sortFileNames(fileNames) 

        # Otherwise we just go through and get the ID's we need
        else:
            for id in idList: 
                fileLoc = '%s%s-log-*%s.dat' % (self.baseDir,self.prefix,id)
                fileNames.extend(glob.glob(fileLoc))

        return fileNames

    # -----------------------------------------------------------------------------
    # def getFileList(self,type,idList=None):
    #     ''' Get a list of input files based on their type, and possibly a number
    #         of unique ID's'''

    #     fileNames = []

    #     # We want all the file names here
    #     if not idList:
    #         lsCommand = 'ls -1 %s%s-%s-%s' % (self.baseDir,self.prefix,type,self.dataName)
    #         fileNames = os.popen(lsCommand).read().split('\n')
    #         fileNames.pop()

    #         # Now sort them
    #         fileNames  = sortFileNames(fileNames) 

    #     # Otherwise we just go through and get the ID's we need
    #     else:
    #         for id in idList: 
    #             lsCommand = 'ls -1 %s%s-log-*%s.dat' % (self.baseDir,self.prefix,id)
    #             fileNames.append(os.popen(lsCommand).read().rstrip('\n'))

    #     return fileNames

# -------------------------------------------------------------------------------
# CLASS SCALAR REDUCE
# -------------------------------------------------------------------------------
class ScalarReduce:
    ''' Helper methods for analzing reduced pimc scalar output data. '''

    # ----------------------------------------------------------------------
    def __init__(self,fileNames,varLabel=None):
        '''Analyze the input files and get all the estimator data.'''

        # Attempt to load numpy
        try:
            import numpy as np
            numpy_loaded = True
        except ImportError:
            numpy_loaded = False 

        # Define the main parameter dictionary
        self.param_ = {}

        # Parameter and estimator name descriptions
        self.descrip = Description()

        # Determine the reduction variable and get its values
        self.reduceLabel = fileNames[0].partition('reduce')[0][-2]
        data = np.loadtxt(fileNames[0],ndmin=2)
        self.param_[self.reduceLabel] = data[:,0]

        # Determine the number of estimators
        self.numEstimators = len(data[0,:])-1

        # Get the estimator column indices
        self.estIndex = getHeadersDict(fileNames[0], removeLab=self.reduceLabel)

        # Extract all other relevant parameters
        for nf,fileName in enumerate(fileNames):
            dataString = fileName.partition('reduce')[-1].rstrip('.dat').split('-')
            while '' in dataString:
                dataString.remove('')
            
            # fill up the param dictionary with values
            for n,dataLabel in enumerate(dataString):
                if dataLabel in self.descrip.paramNames:
                    if nf == 0:
                        self.param_[dataLabel]= [float(dataString[n+1])]
                    else:
                        self.param_[dataLabel].append(float(dataString[n+1]))

        # Now we find the number of unique values of all the parameters
        self.numParams = {}
        for parName,parVals in self.param_.iteritems():
            self.numParams[parName] = len(set(parVals))

        # create an array with the fixed parameters
        self.fixParNames = self.param_.keys()
        self.fixParNames.remove(self.reduceLabel)
        
        # find the name/label of the changing parameter
        if len(fileNames) == 1:
            if varLabel == None:
                self.varLabel = self.fixParNames[0]
            else:
                self.varLabel = varLabel
        else:
            for parName in self.fixParNames:
                if self.numParams[parName] > 1:
                    self.varLabel = parName
                    break;
        
        # Initialize and fill up the main estimator data array
        self.estimator_ = np.zeros([self.numParams[self.varLabel], 
                                   self.numParams[self.reduceLabel], 
                                   self.numEstimators])

        for n,fileName in enumerate(fileNames):
            data = np.loadtxt(fileName,ndmin=2)
            self.estimator_[n,:,: ] = data[:,1:]

    # ----------------------------------------------------------------------
    def getNumVarParams(self):
        '''Return the number of variable parameters.'''
        return self.numParams[self.varLabel]

    # ----------------------------------------------------------------------
    def param(self):
        '''Return the independent parameter over which we are reducing. '''
        return self.param_[self.reduceLabel]

    # ----------------------------------------------------------------------
    def estimator(self,estLabel,ivar):
        '''Return a dependent estimator with a given var number.'''
        return self.estimator_[ivar,:,self.estIndex[estLabel]]

    # ----------------------------------------------------------------------
    def estimatorError(self,estLabel,ivar):
        '''Return a dependent estimator error with a given var number.'''
        return 1.0*self.estimator_[ivar,:,self.estIndex['d_' + estLabel]]

    # ----------------------------------------------------------------------
    def getVarLabel(self,varIndex):
        '''Construct a label for the variable parameter.'''

        labName = self.descrip.paramShortName[self.varLabel]
        labFormat = self.descrip.paramFormat[self.varLabel]
        labValue  = self.param_[self.varLabel][varIndex]
        labUnit = self.descrip.paramUnit[self.varLabel]

        return labName + ' = ' + labFormat % labValue + ' ' + labUnit
#        return lab.rjust(len(lab))

# -------------------------------------------------------------------------------
# CLASS VECTOR REDUCE
# -------------------------------------------------------------------------------
class VectorReduce:
    ''' Helper methods for analzing reduced pimc vector output data. '''

    # ----------------------------------------------------------------------
    def __init__(self,fileNames,estName,varLabel=None):
        '''Analyze the input files and get all the vector estimator data.'''

        # Attempt to load numpy
        try:
            import numpy as np
            numpy_loaded = True
        except ImportError:
            numpy_loaded = False 

        # Define the main parameter dictionary
        self.param_ = {}

        # Parameter and estimator name descriptions
        self.descrip = Description()

        # Determine the reduction variable and get its values
        self.reduceLabel = fileNames[0].partition('reduce')[0][-2]

        # We temporarily load the estimator file to get the values of the reduce
        # variable.  This is easier than globbing it from the vector file
        data = np.loadtxt(fileNames[0].replace(estName,'estimator'),ndmin=2)

        # Get the reduce variable
        self.param_[self.reduceLabel] = data[:,0]

        # Determine the number of estimators
        self.numEstimators = len(data[0,:])-1

        # Extract all other relevant parameters
        for nf,fileName in enumerate(fileNames):
            dataString = fileName.partition('reduce')[-1].rstrip('.dat').split('-')
            while '' in dataString:
                dataString.remove('')
            
            # fill up the param dictionary with values
            for n,dataLabel in enumerate(dataString):
                if dataLabel in self.descrip.paramNames:
                    if nf == 0:
                        self.param_[dataLabel]= [float(dataString[n+1])]
                    else:
                        self.param_[dataLabel].append(float(dataString[n+1]))

        # Now we find the number of unique values of all the parameters
        self.numParams = {}
        for parName,parVals in self.param_.iteritems():
            self.numParams[parName] = len(set(parVals))

        # create an array with the fixed parameters
        self.fixParNames = self.param_.keys()
        self.fixParNames.remove(self.reduceLabel)
        
        # find the name/label of the changing parameter
        if len(fileNames) == 1:
            if varLabel == None:
                self.varLabel = self.fixParNames[0]
            else:
                self.varLabel = varLabel
        else:
            for parName in self.fixParNames:
                if self.numParams[parName] > 1:
                    self.varLabel = parName
                    break;

        # Now we must determine how many vector coordinates there are in our
        # vector estimator
        data = np.loadtxt(fileNames[0],ndmin=2)
        self.numVectorRows = np.size(data,0)
        
        # Initialize and fill up the main estimator data array
        self.estimator_ = np.zeros([self.numParams[self.varLabel], 
                                    self.numVectorRows,
                                    3*self.numParams[self.reduceLabel]]) 

        for n,fileName in enumerate(fileNames):
            data = np.loadtxt(fileName,ndmin=2)
            self.estimator_[n,:,:] = data



    # ----------------------------------------------------------------------
    def getNumVarParams(self):
        '''Return the number of variable parameters.'''
        return self.numParams[self.varLabel]

    # ----------------------------------------------------------------------
    def getNumReduceParams(self):
        '''Return the number of variable parameters.'''
        return self.numParams[self.reduceLabel]

    # ----------------------------------------------------------------------
    def param(self):
        '''Return the independent parameter over which we are reducing. '''
        return self.param_[self.reduceLabel]

#    # ----------------------------------------------------------------------
#    def varParam(self):
#        '''Return the variable parameter array. '''
#        return self.param_[self.reduceLabel]

    # ----------------------------------------------------------------------
    def x(self,varIndex,reduceIndex):
        ''' Return the independent vector variable.

           varIndex:    the index of the variable parameters
           reduceIndex: the index of the reduce parameter
        '''

        return self.estimator_[varIndex,:,3*reduceIndex]

    # ----------------------------------------------------------------------
    def estimator(self,varIndex,reduceIndex):
        ''' Return the dependent vector estimator.

           varIndex:    the index of the variable parameters
           reduceIndex: the index of the reduce parameter
        '''

        return self.estimator_[varIndex,:,3*reduceIndex+1]

    # ----------------------------------------------------------------------
    def estimatorError(self,varIndex,reduceIndex):
        ''' Return the dependent vector estimator error.

           varIndex:    the index of the variable parameters
           reduceIndex: the index of the reduce parameter
        '''

        return self.estimator_[varIndex,:,3*reduceIndex+2]


    # ----------------------------------------------------------------------
    def getReduceLabel(self,reduceIndex):
        '''Construct a label for the reduce parameter.'''

        labName = self.descrip.paramShortName[self.reduceLabel]
        labFormat = self.descrip.paramFormat[self.reduceLabel]
        labValue = self.param_[self.reduceLabel][reduceIndex]
        labUnit = self.descrip.paramUnit[self.reduceLabel]

        return labName + ' = ' + labFormat % labValue + ' ' + labUnit

    # ----------------------------------------------------------------------
    def getVarLabel(self,varIndex):
        '''Construct a label for the variable parameter.'''

        labName = self.descrip.paramShortName[self.varLabel]
        labFormat = self.descrip.paramFormat[self.varLabel]
        labValue  = self.param_[self.varLabel][varIndex]
        labUnit = self.descrip.paramUnit[self.varLabel]

        return labName + ' = ' + labFormat % labValue + ' ' + labUnit


# -------------------------------------------------------------------------------
# CLASS DESCRIPTION
# -------------------------------------------------------------------------------
class Description:
    '''A class which holds descriptions of all the variables used in the path
    ingegral simulations.'''

    # ----------------------------------------------------------------------
    def __init__(self,NDIM=3):
        ''' Defines all maps and dictionaries used in the analysis.'''

        lengthTUnit = r'$[\mathrm{\AA}]$'
        lengthUnit = r'$\mathrm{\AA}$'

        # The name for the density dependent on the dimensionality
        densityName = ['Line', 'Area', 'Volume']

        self.paramNames = ['T','V','u','t','N','n','R','L']

        self.paramShortName = {'T':'T',
                               'V':'V',
                               'u':r'$\mu$',
                               't':r'$\tau$',
                               'N':'N',
                               'n':r'$\rho$',
                               'R':'R',
                               'L':'L'}

        self.paramUnit = {'T':'K',
                          'V':r'$\mathrm{\AA^{%d}}$' % NDIM,
                          'u':'K',
                          't':r'$K^{-1}$',
                          'N':'',
                          'n':r'$\mathrm{\AA}^{-%d}$' % NDIM,
                          'R':lengthUnit,
                          'L':lengthUnit}

        self.paramFormat = {'T':r'%4.2f',
                            'V':r'%3d',
                            'u':r'%+3.1f',
                            't':r'%5.3f',
                            'N':'%3d',
                            'n':r'%f',
                            'R':r'% 4.1f',
                            'L':r'%3d'}

        self.paramLongName = {'T':'Temperature  [K]', 
                              'V':r'Volume  $[\mathrm{\AA^{%d}}]$' % NDIM,
                              'u':'Chemical Potential  [K]', 
                              't':'Imaginary Time Step  [1/K]',
                              'N':'Number of Particles',
                              'n':r'%s Density  $[\mathrm{\AA}^{-%d}]$' % (densityName[NDIM-1],NDIM),
                              'R':'Pore Radius  %s ' % lengthTUnit,
                              'L':'Length %s' % lengthTUnit}

        self.estimatorLongName = {'K':'Kinetic Energy [K]',
                                  'V':'Potential Energy [K]',
                                  'E':'Energy [K]',
                                  'E_mu':r'$E - \mu N$',
                                  'K/N':'Kinetic Energy per Particle [K]',
                                  'V/N':'Potential Energy per Particle [K]',
                                  'E/N':'Energy per Particle [K]',
                                  'N':'Number of Particles',
                                  'N^2':r'(Number of Particles)$^2$',
                                  'density':r'%s Density  $[\mathrm{\AA}^{-%d}]$' % (densityName[NDIM-1],NDIM),
                                  'diagonal':'Diagonal Fraction',
                                  'kappa':r'$\rho^2 \kappa [units]$',
                                  'pair':'Pair Correlation Function [units]',
                                  'radial':r'Radial Density $[\mathrm{\AA}^{-3}]$',
                                  'number':'Number Distribution',
                                  'obdm':'One Body Density Matrix',
                                  'rho_s/rho':'Superfluid Fraction',
                                  'Area_rho_s':'Area Superfluid Fraction',
                                  'rho_s/rho|Z':r'$\rho_s/\rho_0$',
                                  'radwind':r'Radial Winding Superfliud Density $[\mathrm{\AA}^{-3}]$',
                                  'radarea':r'Radial Area Superfliud Density $[\mathrm{\AA}^{-3}]$'
                                 }

        self.estimatorShortName = {'K':'K [K]',
                                  'V':'V [K]',
                                  'E':'E [K]',
                                  'E_mu':r'$E - \mu N$',
                                  'K/N':'K/N [K]',
                                  'V/N':'V/N [K]',
                                  'E/N':'E/N [K]',
                                  'N':'N',
                                  'N^2':r'N$^2$',
                                  'density':r'$\rho [\mathrm{\AA}^{-%d}]$' % (NDIM),
                                  'diagonal':'D',
                                  'kappa':r'$\rho^2 \kappa [units]$',
                                  'pair':'g(r) [units]',
                                  'radial':r'Radial Density $[\mathrm{\AA}^{-3}]$',
                                  'number':'Number Distribution',
                                  'obdm':'One Body Density Matrix',
                                  'rho_s/rho':r'$\rho_s/\rho$',
                                  'rho_s/rho|Z':r'$\rho_s/\rho_0$',
                                  'Area_rho_s':'Area Superfluid Fraction'}

        self.estimatorXLongName = {'number':'Number of Particles',
                                   'pair':'r  %s' % lengthTUnit,
                                   'obdm':'r  %s' % lengthTUnit,
                                   'radial':'r  %s' % lengthTUnit,
                                   'radwind':'r  %s' % lengthTUnit,
                                   'radarea':'r  %s' % lengthTUnit}
