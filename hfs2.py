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
gi = -0.0004
J = array([1/2,  1/2,  3/2])
A = array([2.298157e9, 0.291901e9, 0.050288e9])
Gf = array([1, 1, 1])
gj = array([2, 2/3, 4/3])
baseE = array([0,1e15,1.2e15])
names  = array(['S1/2', 'P1/2', 'P3/2'])
#~ J = array([1/2,  1/2])
#~ A = array([5, 1])
#~ Gf = array([1, 1])
#~ gj = array([2, 2/3])
#~ baseE = array([0,30])
#~ names  = array(['S1/2', 'P1/2'])
B = 1


lev = array([])
print ".... J .... F .... mf .... E .... gF*mF"
for i in range (0,J.size):
    #~ print "--------------"
    #~ print names[i]," hyperfine states"
    F = arange(abs(I-J[i]),I+J[i]+1,1)
    #~ print "F   = ",F
    energy = Ehf(F,I,J[i],A[i])
    #~ print "Ehf = ",energy
    for fi in range(0,F.size):
        mf = arange(-F[fi],F[fi]+1,1)
        #~ print "F  = ",fi	
        #~ print "mf = ",mf	
        for mi in mf:
            
            lev = append(lev,[J[i], F[fi], mi, energy[fi]+baseE[i], landeg(gj[i],gi,F[fi],I,J[i])*mi])


lev = reshape(lev,[-1,5])
print lev


blist = linspace(0,4e8)
ff = zeros((blist.size,lev.shape[0]))
for i in range(0,blist.size):
    for k in range(0,lev.shape[0]):
        ff[i,k] = lev[k,3]+lev[k,4]*blist[i]

for k in range(0,lev.shape[0]):
    plot(blist,ff[:,k])
show()



