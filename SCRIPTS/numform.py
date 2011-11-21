''' Update the number format of output files. '''
import os

# ----------------------------------------------------------------------------
def getPars(fname):
    '''Get the parameters as part of the output file.  
    
    We need to be careful due to the problems associated with splitting at the
    '-' character when there are potential minus signs.
    '''

    dP  = fname.rstrip('.dat').split('-')
    pIndex = []
    for n,d in enumerate(dP):
        if d == '':
            dP[n+1] = '-' + dP[n+1]
            pIndex.append(n)
    for n in pIndex:
        dP.pop(n)

    return dP

# ----------------------------------------------------------------------------
def main():
    outType  = ['estimator','number','obdm','pair','pcycle','super','worm','log','state','radial']

    for type in outType:
        lsCommand = 'ls -1 gce-%s-*' % type
        fileNames = os.popen(lsCommand).read().split('\n')
        fileNames.pop()

        for fname in fileNames:
            s = getPars(fname)
            newName = '%s-%s-%06.3f-%07.3f-%+08.3f-%7.5f-%09d.dat' %\
                    (s[0], s[1], float(s[2]), float(s[3]), float(s[4]),\
                     float(s[5]),int(s[6]))
            os.popen('mv %s %s' % (fname,newName))

if __name__ == "__main__":
    main()
