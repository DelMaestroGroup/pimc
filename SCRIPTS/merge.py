# merge.py
# Adrian Del Maestro
# 10.19.2009
# 
# Merge the results of parallel PIMC output files, and move the 
# originals to an archive

import os,sys,glob
import tarfile
import pimchelp
from optparse import OptionParser
from pylab import *


# -----------------------------------------------------------------------------
def mergeData(pimc,type,newID,skip):
    ''' Merge the results of the PIMC data files to a single file. '''

    fileNames = pimc.getFileList(type)

    numLines = 0
    fileExist = False
    for i,fname in enumerate(fileNames):

        # Does the file exist?
        if len(glob.glob(fname)) > 0:

            # Determine if we are considering an off-diagonal estimator
            diagonalEst = ((fname.find('obdm') == -1) and (fname.find('worm') == -1))

            # Open and prepare the new file
            inFile = open(fname,'r');
            inLines = inFile.readlines();
            inFile.close()
            numLines += (len(inLines)-2) - (diagonalEst)*skip

            if i == 0:
                # get the output file name and open the file for writing
                outName = fname.replace(str(pimc.id[0]),str(newID))
                fileExist = True
                print '%-80s' % outName,
                outFile = open('MERGED/' + outName,'w');

                # replace the old ID for the new one
                inLines[0] = inLines[0].replace(str(pimc.id[0]),str(newID))
                outFile.write(inLines[0])
                outFile.write(inLines[1])

            # strip any comment lines
            inLines.pop(0)
            inLines.pop(0)

            j = 0
            for line in inLines:
                if ((not diagonalEst) or (diagonalEst and j >= skip)):
                    outFile.write(line)
                j += 1

    if fileExist:
        outFile.close() 
        print '%10d' %numLines

    # Now we check if a CYLINDER folder is present, if so, we repeat the process
    if len(glob.glob('CYLINDER')) > 0:
        numLines = 0
        fileExist = False
        for i,fname in enumerate(fileNames):

            cylfname = 'CYLINDER/' + fname
    
            # Does the file exist?
            if len(glob.glob(cylfname)) > 0:
    
                # Determine if we are considering an off-diagonal estimator
                diagonalEst = ((cylfname.find('obdm') == -1) and (cylfname.find('worm') == -1))
    
                # Open and prepare the new file
                inFile = open(cylfname,'r');
                inLines = inFile.readlines();
                inFile.close()
                numLines += (len(inLines)-2) - (diagonalEst)*skip
    
                if i == 0:
                    # get the output file name and open the file for writing
                    outName = cylfname.replace(str(pimc.id[0]),str(newID))
                    fileExist = True
                    print '%-80s' % outName,

                    # We check if we have a CYLINDER directory, if not create it
                    if len(glob.glob('MERGED/CYLINDER')) == 0:
                        os.system('mkdir MERGED/CYLINDER')

                    outFile = open('MERGED/' + outName,'w');
    
                    # replace the old ID for the new one
                    inLines[0] = inLines[0].replace(str(pimc.id[0]),str(newID))
                    outFile.write(inLines[0])
                    outFile.write(inLines[1])
    
                # strip any comment lines
                inLines.pop(0)
                inLines.pop(0)
    
                j = 0
                for line in inLines:
                    if ((not diagonalEst) or (diagonalEst and j >= skip)):
                        outFile.write(line)
                    j += 1
    
        if fileExist:
            outFile.close() 
            print '%10d' %numLines

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = OptionParser() 
    parser.add_option("-T", "--temperature", dest="T", type="float",
                      help="simulation temperature in Kelvin") 
    parser.add_option("-N", "--number-particles", dest="N", type="int",
                      help="number of particles") 
    parser.add_option("-n", "--density", dest="n", type="float",
                      help="number density in Angstroms^{-d}")
    parser.add_option("-t", "--imag-time-step", dest="tau", type="float",
                      help="imaginary time step")
    parser.add_option("-M", "--number-time-slices", dest="M", type="int",
                      help="number of time slices")
    parser.add_option("-u", "--chemical-potential", dest="mu", type="float",
                      help="chemical potential in Kelvin") 
    parser.add_option("-V", "--volume", dest="V", type="float",
                      help="volume in Angstroms^d") 
    parser.add_option("-L", "--Lz", dest="L", type="float",
                      help="Length in Angstroms") 
    parser.add_option("--canonical", action="store_true", dest="canonical",
                      help="are we in the canonical ensemble?")
    parser.add_option("-s", "--skip", dest="skip", type="int",
                      help="how many input lines should we skip?")
    parser.set_defaults(skip=0)

    parser.set_defaults(canonical=False)

    # parse the command line options and get the reduce flag
    (options, args) = parser.parse_args() 
    if len(args) > 0: 
        parser.error("incorrect number of arguments")

    # We check if we have a MERGED directory, if not create it
    if len(glob.glob('MERGED')) == 0:
        os.system('mkdir MERGED')
    
    # Check that we are in the correct ensemble
    pimchelp.checkEnsemble(options.canonical)

    dataName = pimchelp.getFileString(options,reduce=False)

    # Create the PIMC analysis helper and fill up the simulation parameters maps
    pimc = pimchelp.PimcHelp(dataName,options.canonical)
    pimc.getSimulationParameters()

    # We try to find a new PIMCID which is the average of the ones to merge, and
    # make sure it doesn't already exist
    newID = 0
    for id in pimc.id:
        newID += int(id)
    newID = int(newID/(1.0*len(pimc.id)))

    # Now we keep incrementing the ID number until we are sure it is unique
    while ( (len(glob.glob('*estimator*-%09d*' % newID)) > 0) or
           (len(glob.glob('MERGED/*estimator*-%09d*' % newID)) > 0) ):
        newID += 1
    
    # Merge all the output files
    print 'Merged data files:'
    for type in pimc.dataType:
        mergeData(pimc,type,newID,options.skip)

    # copy over the log file
    oldLogName = pimc.getFileList('log')[0]
    newLogName = oldLogName.replace(str(pimc.id[0]),str(newID))
    os.system('cp %s %s' % (oldLogName,'MERGED/'+newLogName))

    # We first create the name of the output tar file
#   mergeName = pimc.getFileList('estimator')[0].rstrip('.dat').replace('estimator','merged')
#   mergeName += '-%09d' % pimc.id[-1]
#   mergeName += '.tar.gz'
#
#   # Archive all the output files  
#   tar = tarfile.open('MERGED/'+mergeName, "w:gz")
#   for type in pimc.outType:
#       fileNames = pimc.getFileList(type)
#
#       # We have to exclude the merged file
#       for i,fname in enumerate(fileNames):
#           if fname.find(str(newID)) != -1:
#               fileNames.pop(i)
#           
#       # Add all the output files to the tar archive
#       for fname in fileNames:
#           tar.add(fname)
#
#           # delete the file that we have now archived
#           os.system('rm %s' % fname)
#
#   tar.close()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

