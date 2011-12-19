# plotoptions.py
# Adrian Del Maestro
# 12.15.2011
# 
# Options that will be used when plotting vector and scalar estimators

import loadgmt,kevent

# -----------------------------------------------------------------------------
def plotOptions(plotType):
    ''' Consctruct a dictionary of plot options based on what type of plot is
    desired.'''

    # First define the common properties
    pOps = {}

    # marker properties
    pOps['markersize'] = 9
    pOps['markeredgewidth'] = 1.0
    pOps['markeredgecolor'] = 'gray'
    pOps['marker'] = 'o'

    # line properties
    pOps['color'] = 'black'
    pOps['linewidth'] = 1.0

    if 'l' in plotType:
        pOps['linestyle'] = '-' 
    else:
        pOps['linestyle'] = 'None' 

    if plotType == 'p':
        pOps['linestyle'] = 'None'

    if plotType == 'l':
        pOps['linewidth'] = 5.0
        pOps['marker'] = None
        pOps['markeredgewidth'] = 0.0
        pOps['markersize'] = 0.0
        pOps.pop('color')

    if 'e' in plotType:
        pOps['capsize'] = 7
        pOps['elinewidth'] = 1.0
        
    return pOps


# -----------------------------------------------------------------------------
def markersColors(numColors):
    '''return the markers and colors used for plotting.'''

    # http://www.graphviz.org/content/color-names
    colors  = loadgmt.getColorList('cb/div','Spectral_05',numColors)
    colors.reverse()
    markers = loadgmt.getMarkerList()

    return markers,colors

