# convert_state.py
# Adrian Del Maestro
# 06.18.2012
# 
# Convert the state file between old and new versions of blitz++ to enable the
# use of legacy input files

import matplotlib

import os,sys
import pyutils
import argparse
import re

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='Convert State Files Between Blitz++ Standards.')
    parser.add_argument('fileNames', help='Scalar estimator files', nargs='+')
    args = parser.parse_args()

    fileNames = args.fileNames

    if len(fileNames) < 1:
        parser.error("Need to specify at least one state file.")

    for fileName in fileNames:
        oldStateFile = file(fileName,'r')
        newState = ""

        foundArray = False
        innerDim = -1
        data = []

        # We first correct for a bug in the way blitz++ outputs arrays which
        # could lead to numbers like %f%f without a space in between them.  This
        # happens for very large or small numbers next to each other.
        oldStateText = oldStateFile.read()
        def splitNumber(matchObj):
            firstNum = matchObj.group(1)
            secondNum = matchObj.group(2)
            if secondNum[0] == '.':
                secondNum = firstNum[-1] + secondNum
                firstNum = firstNum[:-1]
            return  firstNum + ' ' + secondNum
        correctedStateText = re.sub('([-+]?[0-9]*\.[0-9]+)([+-]?[0-9]*\.[0-9]+)',splitNumber,oldStateText)

        # now we test again in case there are three numbers squashed together
        # like %f%f%f
        correctedStateText = re.sub('([-+]?[0-9]*\.[0-9]+)([+-]?[0-9]*\.[0-9]+)',splitNumber,correctedStateText)

        # now we test for any possible numbers in scientific notation
        def splitScientificNumber(matchObj):
            firstNum = matchObj.group('firstNum')
            secondNum = matchObj.group('secondNum')
            if secondNum[0] == '.':
                secondNum = firstNum[-1] + secondNum
                firstNum = firstNum[:-1]
            return  firstNum + ' ' + secondNum

        correctedStateText = re.sub('(?P<firstNum>[-+]?[0-9]*\.[0-9]+([eE][-+]?[0-9][0-9])?)(?P<secondNum>[-+]?[0-9]*\.[0-9]+([eE][-+]?[0-9][0-9])?)',splitScientificNumber,correctedStateText)
        correctedStateText = re.sub('(?P<firstNum>[-+]?[0-9]*\.[0-9]+([eE][-+]?[0-9][0-9])?)(?P<secondNum>[-+]?[0-9]*\.[0-9]+([eE][-+]?[0-9][0-9])?)',splitScientificNumber,correctedStateText)

        # We now iterate through the corrected file
        for line in correctedStateText.splitlines():

            # Skip blank lines or those with a single open/close bracket
            if line.strip() and (line.strip() not in  ['[',']']):

                if foundArray:
                    # Start of array
                    if line.startswith("["):
                        line = line.lstrip("[")

                    # if we find more square opening brackets we have a multi-dimensional array
                    if "[" in line:
                        # determine the dimension of the innermost array
                        line1 = line.lstrip()
                        if innerDim == -1:
                            innerDim = int(line1[0])
                        line1 = line1.replace("%d [" % innerDim,"")
                        line1 = line1.replace("]","")
                        temp = line1.split()
                        data += line1.split()
                    elif line:
                        innerDim = 1
                        line1 = line.lstrip()
                        line1 = line1.replace("]","")
                        data += line1.split()

                    # Create the new data output format
                    if len(data) == innerDim*numRows*numCols:
                        newState += "(0,%d) x (0,%d)\n" % (numRows-1,numCols-1)
                        if innerDim > 1:
                            n = 0
                            for row in range(numRows):
                                if row == 0:
                                    newState += "[ "
                                else:
                                    newState += "  "

                                for col in range(numCols):
                                    newState += "("
                                    for d in range(innerDim):
                                        newState += data[n]
                                        if d < (innerDim-1):
                                            newState += ","
                                        else:
                                            newState += ") "
                                        n += 1
                                if row < numRows-1:
                                    newState += "\n"
                                else:
                                    newState += "]\n\n"
                        else:
                            n = 0
                            for row in range(numRows):
                                if row == 0:
                                    newState += "[ "
                                else:
                                    newState += "  "

                                for col in range(numCols):
                                    newState += "%s " % data[n]
                                    n += 1

                                if row < numRows-1:
                                    newState += "\n"
                                else:
                                    newState += "]\n\n"
                        data = []
                        foundArray = False
                        innerDim = -1

                elif 'x' in line:
                    numRows = int(line.split('x')[0])
                    numCols = int(line.split('x')[1])
                    foundArray = True
                else:
                    newState += line + '\n'

        newState = newState.rstrip()

        # output the new state
        print newState
        oldStateFile.close()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
