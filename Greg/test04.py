#!/usr/bin/env python

from __future__ import division
from numpy import *
from scipy import *
from scipy.linalg import expm
from pylab import *




om = 5
d = 0
G = 6
#~ Listing:  q11, q12, q21, q22
#~ obm = array([[0, -1j*om/2, 1j*om/2, G],\
	                  #~ [0, -1/2*G+1j*d, 0, 1j*om/2],\
	                  #~ [1j*om/2, 0, -1/2*G-1j*d, -1j*om/2],\
	                  #~ [0, 1j*om/2, -1j*om/2, -G]])

obm = array([[0, G, 0, -om],\
	                  [0, -G, 0, om],\
	                  [0, 0, -1/2*G, d],\
	                  [om/2, -om/2, -d, -1/2*G]])


#~ obm[0,:] = array([1,0, 0, 1])
b = array([[1],[0],[0],[0]])

obmin = 1*obm
obmin[0,:] = array([1,1, 0, 0])
bin = array([[1],[0],[0],[0]])
x = linalg.solve(obmin,bin)
#~ print obm,b,x

t = linspace(0,5,1000)
p1 = []
p2 = []
for k in range(0,t.size):
    obmt = expm(obm*t[k])
    p1.append(dot(obmt,b)[0])
    p2.append(dot(obmt,b)[1])
p1 = array(p1)
p2 = array(p2)
plot(t,p1,t,p2,t,p1+p2,'--')
plot(t,p2,[t.min(), t.max()],[x[1], x[1]])
show()


#~ #### Works ####
#~ def pop2(om,d,G):
    #~ obm = array([[0, -1j*om/2, 1j*om/2, G],\
	                  #~ [-1j*om/2, -1/2*G+1j*d, 0, 1j*om/2],\
	                  #~ [1j*om/2, 0, -1/2*G-1j*d, -1j*om/2],\
	                  #~ [0, 1j*om/2, -1j*om/2, -G]])
    #~ obmin = 1*obm
    #~ obmin[0,:] = array([1,0, 0, 1])
    #~ bin = array([[1],[0],[0],[0]])
    #~ x = linalg.solve(obmin,bin)
    #~ return x[3]


#~ d = linspace(-5,5,200)
#~ om = 1
#~ G = 5
#~ p2 = []
#~ for k in range(0,d.shape[0]):
    #~ p2.append(pop2(om,d[k],G))
#~ p2 = array(p2)
#~ plot(d,p2)
#~ show()






