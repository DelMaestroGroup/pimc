''' PimcHelp - helper methods for analzing pimc output data.
'''
import os
from operator import itemgetter, attrgetter

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
    fileTuples = sorted(fileTuples, key=itemgetter(1,2,3,4,5,6))
    
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
    if dataName[0].find('gce') != -1:
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
def getHeadersFromFile(fileName): 
    ''' Get the data column headers from a PIMC output file. '''

    inFile = open(fileName,'r');
    inLines = inFile.readlines();
    if inLines[0].find('PIMCID') != -1:
        headers = inLines[1].split()
    else:
        headers = inLines[0].split()
    headers.pop(0)
    inFile.close()
    return headers

# -------------------------------------------------------------------------------
def checkEnsemble(gce):
    ''' Here we make sure that the correct ensemble flag is specified. '''

    import sys,glob
    gceFiles = (len(glob.glob('gce-*'))>0)
    ceFiles = (len(glob.glob('ce-*'))>0)

    if (gceFiles and not ceFiles) and not gce:
        sys.exit('Need to include -g for the grand canonical ensemble!')

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

    if options.V is not None:
        flagV = "%07.3f" % options.V
        out += '-V-%s' % flagV
    else:
        flagV = "*"

    if options.gce:
        dataName = '%s-%s-%s-%s-*.dat' % (flagT,flagV,flagmu,flagtau)
    else:
        dataName = '%s-%s-%s-%s-*.dat' % (flagT,flagN,flagn,flagtau)

    if reduce:
        outName = '%s-reduce%s' % (options.reduce,out)
        return dataName,outName
    else:
        return dataName

# -------------------------------------------------------------------------------
class PimcHelp:
    ''' Helper methods for analzing pimc output data. '''

    # ----------------------------------------------------------------------
    def __init__(self,dataName,gce,baseDir=''):

        self.dataName = dataName 
        self.baseDir  = baseDir
        self.params = {}
        self.id = []
        if not gce:
            self.prefix='ce'
        else:
            self.prefix='gce'

        # The data file and all output file names
        self.dataType = ['estimator','obdm','pair','pcycle','super','worm','radial']
        if gce:
            self.dataType.append('number')
        self.outType  = ['estimator','number','obdm','pair','pcycle','super','worm','radial','log','state']

    # -----------------------------------------------------------------------------
    def getID(self,fileName): 
        ''' Return the ID number corresponding to a given filename. '''
        #ID = int(fileName.rstrip('.dat').split('-')[-1])
        ID = int(fileName[-13:-4])
        return ID

    # ----------------------------------------------------------------------
    def getParameterMap(self,logName): 
        '''Given a log file name, return the parameter map. '''

        # open the log file
        logFile = open(logName,'r');
        logLines = logFile.readlines();

        # Get the values of all simulation parameters
        paramsMap = {}
        params = False
        for line in logLines:
            if params:
                if line.find('End Simulation Parameters') != -1:
                    params = False
                else:
                    if line.find(':') != -1:
                        keyVal = line.split(':')
                        paramsMap[keyVal[0].strip()] = keyVal[1].strip()
            else:
                if line.find('Begin Simulation Parameters') != -1:
                    params = True

        logFile.close()

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
    def getFileInfo(self,type):
        ''' Get the names of the input files, how many of them there are, and how many
            columns of data they contain.'''

        fileNames = self.getFileList(type)

        # The number of parameter files
        numParams = len(fileNames);

        # Open a sample data file and count the number of columns
        inFile = open(fileNames[0],'r')
        lines = inFile.readlines();
        for line in lines:
            if not (line[0] == '#'):
                cols = line.split()
                break
        numDataCols = len(cols)
        
        return fileNames,numParams,numDataCols

    # -----------------------------------------------------------------------------
    def getFileList(self,type,idList=None):
        ''' Get a list of input files based on their type, and possibly a number
            of unique ID's'''

        fileNames = []

        # We want all the file names here
        if not idList:
            lsCommand = 'ls -1 %s%s-%s-%s' % (self.baseDir,self.prefix,type,self.dataName)
            fileNames = os.popen(lsCommand).read().split('\n')
            fileNames.pop()

            # Now sort them
            fileNames  = sortFileNames(fileNames) 

        # Otherwise we just go through and get the ID's we need
        else:
            for id in idList: 
                lsCommand = 'ls -1 %s%s-log-*%s.dat' % (self.baseDir,self.prefix,id)
                fileNames.append(os.popen(lsCommand).read().rstrip('\n'))

        return fileNames
