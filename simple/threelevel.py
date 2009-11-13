from __future__ import division
import scipy as sp
from scipy.linalg import *

def Mmat(params):
    ''' To evolve the density matrix
    
    Input:
    params : parameters in a standard format =>
      laser1 = (Om1, d1, gg1)
      laser2 = (Om2, d2, gg2)
      params = (G,G1,G2,laser1,laser2)
    
    Output:
    M = density matrix evolving matrix, that is
    dq /dt = M * q
    
    It can be used to calculate complete time evolution.
    E.g. if it is time independent:
    q(t) = q(0) expm(M t) 
    where expm() is the "matrix exponential"
    '''
    (G,G1,G2,laser1,laser2) = params
    (Om1, d1, gg1) = laser1
    (Om2, d2, gg2) = laser2
    # Matrix elements: x11, x22, x33, x12, x13, x23, y12, y13, y23
    M = sp.array([[0, 0, G1, 0, 0, 0, 0, Om1, 0],
                  [0, 0, G2, 0, 0, 0, 0, 0, Om2],
                  [0, 0, -G, 0, 0, 0, 0, -Om1, -Om2],
                  
                  [0, 0, 0, -(gg1+gg2)/2, 0, 0, (d1-d2), Om2/2, Om1/2],
                  [0, 0, 0, 0, -(G1+gg1)/2, 0, Om2/2, d1, 0],
                  [0, 0, 0, 0, 0, -(G2+gg2)/2, -Om1/2, 0, d2],
                                 
                  [0, 0, 0, -(d1-d2), -Om2/2, Om1/2, -(gg1+gg2)/2, 0, 0],
                  [-Om1/2, 0, Om1/2, -Om2/2, -d1, 0, 0, -(G1+gg1)/2, 0],
                  [0, -Om2/2, Om2/2, -Om1/2, 0, -d2, 0, 0, -(G2+gg2)/2]])
    return M

def timeseries(tseries, params, q0):
    ''' Time-dependent state evolution
        
    Input parameters:
    tseries = (tstart, tend, tnum) : time points to check (units of 1/G)
    params = parameters in the standard format
    q0 = starting configuration
    
    Outout:
    (t, p1, p2, p3) : vectors of times and the three populations

    Example:
    laser1 = (Om1, d1, gg1)
    laser2 = (Om2, d2, gg2)
    params = (G,G1,G2,laser1,laser2)
    (t, p1, p2, p3) = timeseries((0,20,101),params, q0)
    # This will calculate the time evolution between 0 and 20/G,
    # using starting configuration q0
    '''
    # Setting up parameters
    (tstart, tend, tnum) = tseries
    t = sp.linspace(tstart, tend , tnum)
    p1 = sp.zeros((len(t),1))
    p2 = sp.zeros((len(t),1))
    p3 = sp.zeros((len(t),1))
    # Time evolution
    Mnow = Mmat(params)
    for i in range(0, len(t)):
        temp = sp.dot(expm(Mnow * t[i]), q0)
        p1[i] = temp[0]
        p2[i] = temp[1]
        p3[i] = temp[2]
    return (t, p1, p2, p3)

def spectra(detseries, laser, params, q0):
    ''' Spectra calculation
    
    Input parameters:
    detseries = (detstart, detend, detnum) : detunings to check (units of G)
    laser = 1/2 : which laser to scan
    params = parameters in the standard format
    q0 = starting configuration
    
    Outout:
    (det, p1, p2, p3) : vectors of detunings and the three populations

    Example:
    laser1 = (Om1, d1, gg1)
    laser2 = (Om2, d2, gg2)
    params = (G,G1,G2,laser1,laser2)
    (det, p1, p2, p3) = spectra((-10,10,401), 1, params, q0)
    # This will calculate the spectra between -10G and 10G detuning,
    # scanning laser1, using starting configuration q0
    '''
    # Setting up parameters
    (detstart, detend, detnum) = detseries
    det = sp.linspace(detstart, detend , detnum)
    p1 = sp.zeros((len(det),1))
    p2 = sp.zeros((len(det),1))
    p3 = sp.zeros((len(det),1))
    (G,G1,G2,laser1,laser2) = params
    (Om1, d1, gg1) = laser1
    (Om2, d2, gg2) = laser2
    # A very long time to go ahead and check the "final" population
    longt = 1000
    # Scan the laser detuning
    for i in range(0, len(det)):
        if laser == 1 :
            laser1 = (Om1, det[i], gg1)
        else :
            laser2 = (Om2, det[i], gg2)
        params = (G,G1,G2,laser1,laser2)
        temp = sp.dot(expm(Mmat(params) * longt), q0)
        p1[i] = temp[0]
        p2[i] = temp[1]
        p3[i] = temp[2]
    return (det, p1, p2, p3)
