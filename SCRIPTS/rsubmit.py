#!/usr/bin/python 
#
# rsubmit.py
# Adrian Del Maestro
# 09.04.2009
#
# Generate a torque submit script for the pimc code which creates
# a pbs file for various sets of parameters.  This restart script
# reads all log files in a current directory by default, or will
# take a list of id numbers from the command line

import sys,os
import pimchelp
from optparse import OptionParser

# -----------------------------------------------------------------------------
def getPIMCommand(fname):
    ''' Get the command string from the log file.'''
    # Open the log file and get the command string
    logFile = open(fname,'r')
    logLines = logFile.readlines()
    line = logLines[2]

    return line[2:]

# -----------------------------------------------------------------------------
def westgrid(logFileNames,outName):
    ''' Write a pbs submit script for westgrid. '''

    # Open the pbs file and write its header
    fileName = 'resubmit-pimc%s.pbs' % outName
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash\n
#PBS -l walltime=120:00:00
#PBS -N PIMC
#PBS -V
#PBS -j oe
#PBS -o out/pimc-${PBS_JOBID}\n
# Do not send email
#PBS -M agdelma@phas.ubc.ca
#PBS -m n\n
# Start job script
cd $PBS_O_WORKDIR
echo \"Starting run at: `date`\"\n''')
    
    numID = len(logFileNames)

    if numID > 1:
        pbsFile.write('''\ncase ${PBS_ARRAYID} in\n''')
    
    # Create the command string and make the case structure
    for n,fname in enumerate(logFileNames):
        # Get the command string
        command = getPIMCommand(fname)
        if numID > 1:
            pbsFile.write('%d)\n%s;;\n' % (n,command))
        else:
            pbsFile.write(command)
    
    if numID > 1:
        pbsFile.write('esac\necho \"Finished run at: `date`\"')
        print '\nSubmit jobs with: qsub -t 0-%d %s\n' % (numID-1,fileName)
    else:
        pbsFile.write('echo \"Finished run at: `date`\"')
        print '\nSubmit jobs with: qsub %s\n' % fileName
    pbsFile.close();

# -----------------------------------------------------------------------------
def sharcnet(logFileNames,outName):
    ''' Write a pbs submit script for sharcnet. '''

    # Open the script file and write its header
    fileName = 'resubmit-pimc%s' % outName
    scriptFile = open(fileName,'w')
    scriptFile.write('''#!/bin/bash
# Sharcnet pimc submit script\n\n''')

    # Get the command string and output to submit file
    name = '''out/pimc-%J'''
    for n,fname in enumerate(logFileNames):
        command = getPIMCommand(fname)
        scriptFile.write('sqsub -q serial -o %s -r 6d %s' % (name,command))
    scriptFile.close();
    os.system('chmod u+x %s'%fileName)

# -----------------------------------------------------------------------------
def scinet(logFileNames,outName):
    ''' Write a pbs submit script for scinet. '''

    # Open the pbs file and write its header
    fileName = 'resubmit-pimc%s.pbs' % outName
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash
# MOAB/Torque submission script for multiple serial jobs on SciNet GPC
#
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N serialx8_pimc
#PBS -V
#PBS -j oe
#PBS -o out/pimc-${PBS_JOBID}
# Do not send email
#PBS -M agdelma@phas.ubc.ca
#PBS -m n

# Start job script
cd $PBS_O_WORKDIR
echo \"Starting run at: `date`\"\n\n''')

    # Get the command string and output to submit file
    for n,fname in enumerate(logFileNames):
        command = getPIMCommand(fname).rstrip('\n')
        pbsFile.write('(%s) &\n' % command)
    pbsFile.write('wait')
    pbsFile.close();
    print '\nSubmit job with: qsub %s\n' % fileName

# -----------------------------------------------------------------------------
# Begin Main Program
# -----------------------------------------------------------------------------
def main():

    # setup the command line parser options 
    parser = OptionParser() 
    parser.add_option("-T", "--temperature", dest="T", type="float", \
            help="simulation temperature in Kelvin") 
    parser.add_option("-N", "--number-particles", dest="N", type="int",\
            help="number of particles") 
    parser.add_option("-n", "--density", dest="n", type="float",\
            help="number density in Angstroms^{-d}")
    parser.add_option("-P", "--number-time-slices", dest="P", type="int",\
            help="number of time slices")
    parser.add_option("-u", "--chemical-potential", dest="mu", type="float",\
            help="chemical potential in Kelvin") 
    parser.add_option("-L", "--length", dest="L", type="float",\
            help="length in Angstroms") 
    parser.add_option("-t", "--imag-time-step", dest="tau", type="float",\
            help="imaginary time step") 
    parser.add_option("--canonical", action="store_true", dest="canonical", 
                      help="are we in the canonical ensemble?")
    parser.add_option("-i", "--id", action="append", dest="pimcID", type="int",\
            help="a list of PIMC ID numbers to include")
    parser.add_option("-e", "--exclude", action="append", dest="exID", type="int",\
            help="a list of PIMC ID numbers to exclude")
    parser.add_option("-c", "--cluster", dest="cluster", choices=['westgrid','sharcnet','scinet'],\
            help="target cluster: [westgrid,sharcnet,scinet]") 

    parser.set_defaults(canonical=False)

    # parse the command line options and get the reduce flag
    (options, args) = parser.parse_args() 

    if (not options.cluster):
        parser.error("need to specify a cluster")

    # Check that we are in the correct ensemble
    pimchelp.checkEnsemble(options.canonical)

    # create a file string that will be used to name the submit file
    outName = ''
    if options.T:
        outName += '-%06.3f' % options.T
    if options.L:
        outName += '-%07.3f' % options.L

    # Get the data string and create the pimc helper object
    dataName = pimchelp.getFileString(options,reduce=False)
    pimc = pimchelp.PimcHelp(dataName,options.canonical)

    # We get either all the log files in the current directory, or just the
    # requested files by their ID number
    logFileNames  = pimc.getFileList('log',idList=options.pimcID)

    # If we have excluded any ID's we remove them from the list
    if options.exID:
        for id in options.exID:
            for n,fname in enumerate(logFileNames):
                if int(id) == pimc.getID(fname):
                    logFileNames.pop(n)

    # Now create the submission files
    if options.cluster == 'westgrid':
        westgrid(logFileNames,outName)

    if options.cluster == 'sharcnet':
        sharcnet(logFileNames,outName)

    if options.cluster == 'scinet':
        scinet(logFileNames,outName)

# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
