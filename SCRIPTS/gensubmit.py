#!/usr/bin/python 
#
# gensubmit.py
# Adrian Del Maestro
# 06.03.2009
#
# Generate a torque submit script for the pimc code which creates
# a pbs file for various sets of parameters.

import os,sys
from optparse import OptionParser

# -----------------------------------------------------------------------------
def westgrid(staticPIMCOps,numOptions,optionValue,outName):
    ''' Write a pbs submit script for westgrid. '''

    # Open the pbs file and write its header
    fileName = 'submit-pimc%s.pbs' % outName
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
echo \"Starting run at: `date`\"

case ${PBS_ARRAYID} in\n''')

    # Create the command string and make the case structure
    for n in range(numOptions):
        if (optionValue.has_key('p') or staticPIMCOps.find('-p') != -1):
            command = './pimc.e '
        else:
            command = './pimc.e -p %d ' % (n)

        for flag,val in optionValue.iteritems():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps
        pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,2*n,command))
    
    pbsFile.write('esac\necho \"Finished run at: `date`\"')
    pbsFile.close();
    
    print '\nSubmit jobs with: qsub -t 0-%d %s\n' % (numOptions-1,fileName)

# -----------------------------------------------------------------------------
def sharcnet(staticPIMCOps,numOptions,optionValue,outName):
    ''' Write a submit script for sharcnet. '''

    # Open the script file and write its header
    fileName = 'submit-pimc%s' % outName
    scriptFile = open(fileName,'w')
    scriptFile.write('''#!/bin/bash
# Sharcnet pimc submit script\n\n''')

    # Create the command string and output to submit file
    for n in range(numOptions):
        name = '''out/pimc-%J'''
        if optionValue.has_key('p'):
            command = './pimc.e '
        else:
            command = './pimc.e -p %d ' % (n)
            
        for flag,val in optionValue.iteritems():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps
        scriptFile.write('sleep %d\nsqsub -q serial -o %s --mpp=1G -r 6d %s\n' % (10,name,command))
    scriptFile.close();
    os.system('chmod u+x %s'%fileName)

# -----------------------------------------------------------------------------
def scinet(staticPIMCOps,numOptions,optionValue,outName):
    ''' Write a pbs submit script for scinet. '''

    if numOptions != 8:
        print 'For scinet, must submit in multiples of 8'
        return

    # Open the pbs file and write its header
    fileName = 'submit-pimc%s.pbs' % outName
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

    # Create the command string and output to submit file
    for n in range(numOptions):
        command = './pimc.e '
        for flag,val in optionValue.iteritems():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps.rstrip('\n')
        pbsFile.write('(sleep %02d; %s) &\n' % (10,command))
    pbsFile.write('wait')
    pbsFile.close();
    print '\nSubmit job with: qsub %s\n'%fileName

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = OptionParser() 
    parser.add_option("-c", "--cluster", dest="cluster", choices=['westgrid','sharcnet','scinet'],\
            help="target cluster: [westgrid,sharcnet,scinet]") 

    # parse the command line options
    (options, args) = parser.parse_args() 
    inFileName = args[0]

    if (not options.cluster):
        parser.error("need to specify a cluster")

    # We open up the input file, and read in all lines.
    inFile = open(inFileName,'r')
    inLines = inFile.readlines();

    # The first line of the file contains all static pimc options
    staticPIMCOps = inLines[0].rstrip('\n')

    # The next lines contains the short-name of the pimc option
    # and a list of values.  
    optionValue = {}
    numOptions = []
    for line in inLines[1:]:
        option = line.split()
        flag = option.pop(0)
        # Now we determine if we have a range or if we have actually included the values
        if option[0].find(':') != -1:

            # determine if we have floats or ints
            if option[0].find('.') != -1:
                type_converter = lambda x: float(x)
            else:
                type_converter = lambda x: int(x)

            dataRange = option[0].split(':')

            # Parse the range string
            dataMin = type_converter(dataRange[0])
            dataMax = type_converter(dataRange[1])
            dataStep = type_converter(dataRange[2])

            # construct the option list
            vals = [dataMin]
            while vals[-1] < dataMax:
                vals.append(vals[-1]+dataStep)

            # assign the option list
            optionValue[flag] = vals
            numOptions.append(len(vals))
        else:
            # store the typed out values
            optionValue[flag] = option
            numOptions.append(len(option))
    
    # We try to extract the temperature and volume to make an output string
    outName = ''
    findInput = staticPIMCOps.split();
    n = 0
    for input in findInput:
        if input == '-T':
            outName += '-%06.3f' % float(findInput[n+1])
            break
        n += 1
    n = 0
    for input in findInput:
        if input == '-L':
            outName += '-%07.3f' % float(findInput[n+1])
            break
        n += 1

    # We make sure all option strings have the same length 
    if len(numOptions) > 0:
        for n in numOptions:
            if n != numOptions[0]:
                print 'Not all parameters have the same number of values!'
                sys.exit()
        numOptions = numOptions[0]

    if options.cluster == 'westgrid':
        westgrid(staticPIMCOps,numOptions,optionValue,outName)

    if options.cluster == 'sharcnet':
        sharcnet(staticPIMCOps,numOptions,optionValue,outName)

    if options.cluster == 'scinet':
        scinet(staticPIMCOps,numOptions,optionValue,outName)

# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
