#!/usr/bin/env python
from __future__ import division
from numpy import *
from pylab import *
from scipy import *
from physicspy.quantum import *

def Ehf(F,I,J,A):
	return A/2*(F*(F+1) - I*(I+1) - J*(J+1))

def gf(F,J,I):
    return (F*(F+1) + J*(J+1) - I*(I+1))/(2*F*(F+1))




I = 7/2
gI = -0.00039885395
gS = 2.0023193043622
gL = 0.99999587
uB = 1.399624e6
J = array([1/2, 1/2, 3/2])
#~ J = array([1/2])
Ahfs = array([2.2981572425e9,  291.9201e6, 50.28827e6])
Bhfs = array([0, 0, -0.4934e6])
Chfs = array([0, 0, 0.560e3])
gJ = array([2.00254032, 0.665900, 1.334])
baseE = array([0, 0, 0])
names  = array(['S1/2', 'P1/2', 'P3/2'])


lev = array([])
print "J....F....mf....dEhfs(MHz)....dE(MHz)/G"
for i in range (0,J.size):
    #~ print "--------------"
    #~ print names[i]," hyperfine states"
    F = arange(abs(I-J[i]),I+J[i]+1,1)
    #~ print "F   = ",F
    #~ energy = Ehf(F,I,J[i],Ahfs[i])
    #~ print "Ehf = ",energy
    for fi in range(0,F.size):
        energy = energyhfsalkali(F[fi],J[i],I,Ahfs[i],Bhfs[i],Chfs[i])
        gf = landeg(gJ[i],gI,F[fi],I,J[i])
        EB = uB*gf
        mf = arange(-F[fi],F[fi]+1,1)
        for mi in mf:
            lev = append(lev,[J[i], F[fi], mi, energy+baseE[i], EB*mi])


lev = reshape(lev,[-1,5])
lev[:,3] /= 1e6
lev[:,4] /= 1e6
print lev[:,2:]
