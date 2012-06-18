#!/usr/bin/python
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
        for line in oldStateFile.readlines():

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
                    line1 = line1.replace("-"," -")
                    temp = line1.split()

                    # Correct for poor blitz output which could allow two numbers to
                    # be included in single record
                    numList = []
                    for n,t in enumerate(temp):
                        if t.count('.') == 2:
                            decPosition = t[::-1].index('.')
                            secondNum = t[decPosition-1:]
                            firstNum = t.rstrip(secondNum)
                            numList.append(firstNum)
                            numList.append(secondNum)
                        else:
                            numList.append(t)

                    data += numList
                else:
                    innerDim = 1
                    line1 = line.lstrip()
                    line1 = line1.replace("]","")
                    data += line1.split()

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
                                newState += "]\n"
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
                                newState += "]\n"
                    data = []
                    foundArray = False
                    innerDim = -1

            elif 'x' in line:
                numRows = int(line.split('x')[0])
                numCols = int(line.split('x')[1])
                foundArray = True
            else:
                newState += line

        newState = newState.rstrip()
        print newState
        oldStateFile.close()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
