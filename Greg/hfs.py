#!/usr/bin/env python
from __future__ import division
from numpy import *
from pylab import *

def Ehf(F,I,J,A):
	return A/2*(F*(F+1) - I*(I+1) - J*(J+1))

I = 7/2
J = array([1/2,  1/2,  3/2])
A = array([1, 1, 1])
Gf = array([1, 1, 1])
names  = array(['S1/2', 'P1/2', 'P3/2'])
B = 1

for i in range (0,J.size):
    print "--------------"
    print names[i]," hyperfine states"
    F = arange(abs(I-J[i]),I+J[i]+1,1)
    print "F   = ",F
    energy = Ehf(F,I,J[i],A[i])
    print "Ehf = ",energy
    for fi in F:
        mf = arange(-fi,fi+1,1)
        print "F  = ",fi	
        print "mf = ",mf	

