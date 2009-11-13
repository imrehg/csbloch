from __future__ import division
import scipy as sp
from scipy.linalg import *

# Helper functions
def conn(j, k, n):
    ''' Iterate through the non-diagonal elements '''
    if (j == k) :
        return -1
    if (j > k):
        j, k = k, j
    out = -1
    for temp in range(0, j+1):
        out += (n - 1) - temp
    out -= (n - 1) - k
    return out

def xv(j, k, n):
    ''' Find x(j, k) serial number in list:
    Order: first non diagonal elements (0,0)...(n,n) then off-diagonals
    '''
    if (j >= n) | (k >= n):
        return -1
    if (j == k):
        # Diagonal elements
        return j
    else :
        # Coherence, real part
        out = conn(j, k, n) + n
    return out

def yv(j, k, n):
    ''' Find y(j, k) serial number in list:
    Order: Off-diagonals only, after x values.
    
    Return:
    (y, mul) :  y = -1 if invalid (j, k) values
                mul = -1 if (j > k), otherwise =1
    '''
    if (j >= n) | (k >= n):
        return (-1, -1)
    if (j == k):
        # Diagonal elements
        return (-1, -1)
    else :
        # Coherence, imaginary part
        multi = 1
        if (j > k):
            multi = -1
        out = conn(j, k, n) + 2 * n
    return (out, multi)

def ChooseMat(j, k, n, Choose):
    ''' Choose an element from a matrix in off-diagonal elements order
    Useful to handle state-state interactions
    '''
    iter = conn(j, k, n)
    if (iter < 0):
        return 0
    else:
        return Choose[iter]

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
    # States:
    n = 3
    # ChooseMat: offdiagonals
    Om = [0, Om1/2, Om2/2]
    A = [0, G1, G2]
    Detu = [d1-d2, d1, d2]
    # Diagonals
    # Gamma - total decay
    Gam = [0, 0, G]
    # Big gamma, width of states
    GM = sp.array([0, 0, G])
    # Laser linewidth
    lw = [gg1+gg2, gg1, gg2]
    
    # Matrix elements: x11, x22, x33, x12, x13, x23, y12, y13, y23
    Mm = sp.zeros((n**2, n**2))

    # Diagonal elements
    for j in range(0, n):
        for l in range(0, 3):
            y, mul = yv(j, l, n)
            if (y > -1):
                M = ChooseMat(l, j, n, Om)
                Mm[j, y] += 2 * M * mul
            if (j < l):
                Mm[j, l] += ChooseMat(j, l, n, A)
        Mm[j, j] += -Gam[j]

    # Off-diagonal, real part 
    for j in range(0, n):
        for k in range(j+1,n):
            x = xv(j, k, n)
            for l in range(0, n):
                y, mul = yv(l, k, n)
                if (y > -1):
                    M = ChooseMat(j, l, n, Om)
                    Mm[x, y] += -M * mul
                y, mul = yv(j, l, n)
                if (y > -1):
                    M = ChooseMat(l, k, n, Om)
                    Mm[x, y] += M * mul
            Mm[x, x] += -(GM[j] + GM[k] + ChooseMat(j, k, n, lw)) / 2
            y, mul = yv(j, k, n)
            if (y > -1):
                Mm[x, y] += mul * ChooseMat(j, k, n, Detu)

    # Off-diagonal, imaginary part 
    for j in range(0, n):
        for k in range(j+1,n):
            y, mul = yv(j, k, n)
            for l in range(0, n):
                x = xv(l, k, n)
                M = ChooseMat(j, l, n, Om)
                Mm[y, x] += M
                x = xv(j, l, n)
                M = ChooseMat(l, k, n, Om)
                Mm[y, x] += -M
            Mm[y, y] += -(GM[j] + GM[k] + ChooseMat(j, k, n, lw)) / 2
            x = xv(j, k, n)
            Mm[y, x] += - ChooseMat(j, k, n, Detu)
    return Mm

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
