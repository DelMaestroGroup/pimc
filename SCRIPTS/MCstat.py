#!/usr/bin/python
# pimcstat.py
# Chris Herdman
# 11.30.2012
# 
# Statistical analysis library for Monte Carlo Data
import math as m
import numpy as np
import os.path
import sys


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Begin bin function
# -----------------------------------------------------------------------------
# Returns error as a funciton of bin level (D)
# Requires raw Monte Carlo bin data in horizontal rows
# -----------------------------------------------------------------------------
def bin(MC):
    
    # minimum number of MC bins required
    min_bin = 8
    # initialize B to MC data
    B = MC
    # Resize if 1D array
    if B.ndim == 1:
        B.resize(B.shape[0],1)
    # Define number of binning levels
    Nl = int(m.floor(m.log(B.shape[0]/min_bin,2))+1)
    # initialize D
    D = np.zeros((Nl,B.shape[1]))
    
    # First level of binning is raw data
    D[0,:] = np.std(MC,0)/m.sqrt(B.shape[0]-1)
    
    # Binning loop over levels l
    for l in range(1,Nl):
        
        # Bin pairs of bins: if odd # of bins, truncate first bin
        if ((B.shape[0] % 2) == 0):
            B = (B[::2,:]+ B[1::2,:])/2
        else:
            B = (B[1::2,:]+ B[2::2,:])/2
        
        # Error at level l
        D[l,:] = np.std(B,0)/m.sqrt(B.shape[0]-1)
    
    return D
# -----------------------------------------------------------------------------
# end of bin function
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Begin bin_conv function
# -----------------------------------------------------------------------------
# Returns convergence factor (CF) of binning analysis
# Requires processed bins from bin function
# -----------------------------------------------------------------------------
def bin_conv(Delta,dDelta=None):
    
    if Delta.ndim == 1:
        Delta.resize((len(Delta),1))
        if dDelta is not None:
            dDelta.resize((len(dDelta),1))
    
    # initialize convergence factor
    CF = np.zeros(Delta.shape[1])
    dCF = np.zeros(Delta.shape[1])
    
    #Loop of columns of D to compute tau,CF
    for c in range(Delta.shape[1]):
        if Delta[-2,c] != 0.0:
            CF[c] = 1.0-(Delta[-1,c]-Delta[-2,c])/Delta[-2,c]
            if dDelta is not None:
                dCF[c] = np.sqrt( ((1-CF[c])**2)*( \
                                    (dDelta[-1,c]**2+dDelta[-2,c]**2)/\
                                    ((Delta[-1,c]-Delta[-2,c])**2) +\
                                    (dDelta[-2,c]**2)/(Delta[-2,c]**2) ) )
    
    return {'CF':CF,'dCF':dCF}
# -----------------------------------------------------------------------------
# end of bin_conv function
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Begin bin_ac function
# -----------------------------------------------------------------------------
# Returns estimate of autocorrelation time in bins from binning analysis (tau)
# Requires processed bins from bin function
# -----------------------------------------------------------------------------
def bin_ac(D,dD=None):
    
    if D.ndim == 1:
        D.resize((len(D),1))
        if dD is not None:
            dD.resize((len(dD),1))
        
    # initialize autocorrelation time
    tau = np.zeros(D.shape[1])
    dtau = np.zeros(D.shape[1])
    
    #Loop of columns of D to compute tau,CF
    for c in range(D.shape[1]):
        if D[0,c] != 0.0:
            tau[c] = 0.5*((D[-1,c]/D[0,c])**2-1)
            if dD is not None:
                dtau[c] = np.sqrt(tau[c]**2*( (dD[-1,c]/D[-1,c])**2 + \
                                   (dD[0,c]/D[0,c])**2 ) )
            
    return {'tau':tau,'dtau':dtau}
# -----------------------------------------------------------------------------
# end of bin_ac function
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Begin autocorr function
# -----------------------------------------------------------------------------
# Returns autocorrelation as a funciton of MC update time
# Averages dt over entire MC run, from dt=1 to dt=T/2
# -----------------------------------------------------------------------------
def autocorr(MC):
    
    # Reshape MC to 2D array if column vector
    if MC.ndim == 1:
        MC.resize(MC.shape[0],1)
    
    # Number of MC bins
    T = MC.shape[0]
    
    # Number of AC steps
    Nt = T/2+1
    
    ac = np.zeros((Nt,MC.shape[1]))
    dac = np.zeros((Nt,MC.shape[1]))
    
    A2 = np.average(MC,0)**2
    dA2 = np.sqrt(2*A2)*np.std(MC,0)/np.sqrt(T)
    AA = np.average(MC*MC,0)
    dAA = np.std(MC*MC,0)/np.sqrt(T)
    ac[0,:] = np.ones((1,MC.shape[1]))
    for dt in range(1,Nt):
        AtA0 = MC[:-dt,:]*MC[dt:,:]
        ac[dt,:] = (np.average(AtA0,0)-A2)/(AA-A2)
        dac[dt,:] = np.sqrt( (ac[dt,:]**2)*( \
                    (np.std(AtA0,0)**2/AtA0.shape[0]+dA2**2)/ \
                    (np.average(AtA0,0)-A2)**2 + (dAA**2+dA2**2)/(AA-A2)**2) )
    
    return {'correlation':ac, 'error':dac}
    
# -----------------------------------------------------------------------------
# end of autocorr function
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
