#!/usr/bin/env python

from __future__ import division
from numpy import *
from scipy import *
from scipy.linalg import expm
from pylab import *




om = 0.1
d = 1
G = 1.5
obm = array([[0, G, -1j*om/2, 1j*om/2],\
	    [1j*om/2, -1j*om/2, -1j*d, 0],\
	    [-1j*om/2, 1j*om/2, 0, 1j*d],\
	    [0, -G, 1j*om/2, 1j*om/2]])
#~ obm[0,:] = array([1, 1, 0, 0])
b = array([[1],[0],[0],[0]])
#~ x = linalg.solve(obm,b)
#~ print obm,b,x
t = 0

t = linspace(0,10)
p1 = []
p2 = []
for k in range(0,t.size):
    obmt = expm(obm*t[k])
    p1.append(dot(obmt,b)[0])
    p2.append(dot(obmt,b)[1])
p1 = array(p1)
p2 = array(p2)
plot(t,p1,t,p2,t,p1+p2,'--')
show()