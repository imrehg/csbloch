from __future__ import division
import scipy as sp
from scipy.linalg import *
import pylab as pl

##Level structure
##
##     3
##    ---
##
##
##       ---
## ---    2
##  1

# Excited state lifetimes
# Total should be 1, everything else is scaled with it
G = 1.0
G1 = 0.5
G2 = G - G1
# Rabi frequencies
Om1 = 1.0
Om2 = 0.5
# Laser linewidths
gg1 = 0
gg2 = gg1
# Detuning
d1 = 0
d2 = 0
# Starting density matrix elements
q0 = sp.zeros((9,1))
q0[0] = 1

def Mmat(params):
    ''' Liuville matrix for the density matrix '''
    (G,G1,G2,laser1,laser2) = params
    (Om1, d1, gg1) = laser1
    (Om2, d2, gg2) = laser2
    ##Matrix elements: p1, p2, p3, q12, q21, q13, q31, q23, q32
    M = sp.array([[0, 0, G1, 0, 0, 1j*(-Om1/2), 1j*(Om1/2), 0, 0],
                  [0, 0, G2, 0, 0, 0, 0, 1j*(-Om2/2), 1j*(Om2/2)],
                  [0, 0, -G, 0, 0, 1j*(Om1/2), 1j*(-Om1/2), 1j*(Om2/2), 1j*(-Om2/2)],
                  [0, 0, 0, 1j*(d1-d2), 0, 1j*(-Om2/2), 0, 0, 1j*(Om1/2)],
                  [0, 0, 0, 0, 1j*(-d1+d2), 0, 1j*(Om2/2), 1j*(-Om1/2), 0],
                  [1j*(-Om1/2), 0, 1j*(Om1/2), 1j*(-Om2/2), 0, -(G1+gg1)/2-1j*d1, 0, 0, 0],
                  [1j*(Om1/2), 0, 1j*(-Om1/2), 0, 1j*(Om2/2), 0, -(G1+gg1)/2+1j*d1, 0, 0],
                  [0, 1j*(-Om2/2), 1j*(Om2/2), 0, 1j*(-Om1/2), 0, 0, -(G2+gg2)/2-1j*d2, 0],
                  [0, 1j*(Om2/2), 1j*(-Om2/2), 1j*(Om1/2), 0, 0, 0, 0, -(G2+gg2)/2+1j*d2]])
    return M

def timeseries(tseries, params, q0):
    ''' Time-dependent state evolution '''
    (tstart, tend, tnum) = tseries
    t = sp.linspace(tstart, tend , tnum)
    p1 = sp.zeros((len(t),1))
    p2 = sp.zeros((len(t),1))
    p3 = sp.zeros((len(t),1))   
    for i in range(0, len(t)):
        temp = sp.dot(expm(Mmat(params) * t[i]), q0)
        p1[i] = temp[0]
        p2[i] = temp[1]
        p3[i] = temp[2]
    return (t, p1, p2, p3)
    
laser1 = (Om1, d1, gg1)
laser2 = (Om2, d2, gg2)
params = (G,G1,G2,laser1,laser2)

(t, p1, p2, p3) = timeseries((0,10,101),params, q0)

pl.plot(t, p1, 'b--', label='Ground state 1')
pl.plot(t, p2, 'g-.', label='Ground state 2')
pl.plot(t, p3, 'r-', label='Excited state 3')
pl.xlabel('time')
pl.ylabel('Population')
pl.ylim([0,1])
pl.legend(loc="best")
pl.show()
