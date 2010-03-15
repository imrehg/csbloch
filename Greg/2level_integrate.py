#!/usr/bin/env python

from __future__ import division
from numpy import *
from scipy import *
from scipy.linalg import expm
from pylab import *

om = 10
d = 0
G = 0.1

obm = array([[0, G, 0, -om], \
             [0, -G, 0, om], \
             [0, 0, -1/2*G, d], \
             [om/2, -om/2, -d, -1/2*G]])


b = array([[1],[0],[0],[0]])

obmin = 1*obm
obmin[0,:] = array([1,1, 0, 0])
bin = array([[1],[0],[0],[0]])
x = linalg.solve(obmin,bin)

t = linspace(0,5,1000)
p1 = []
p2 = []
for k in range(0,t.size):
    obmt = expm(obm*t[k])
    p1.append(dot(obmt,b)[0])
    p2.append(dot(obmt,b)[1])
p1 = array(p1)
p2 = array(p2)





### Integrate directly

from scipy.integrate import ode

def f(t0, y):
    y1 = zeros(4)
    for i in range(0, 4):
        y1[i] += dot(obm[i,:],y)
    return(y1)

y0 = b
t0 = 0
dt = 1
tend = 5
r = ode(f).set_integrator('zvode', with_jacobian=False)
r.set_initial_value(y0, t0)
rt = array([t0])
ry = array(y0.T)
while r.successful() and r.t < tend:
    r.integrate(r.t+dt)
    #print ry
    #print r.y
    #print r.t
    rt = append(rt, r.t)
    ry = append(ry, r.y.T, 0)

plot(t,p1,t,p2,t,p1+p2,'--')
plot(t,p2,[t.min(), t.max()],[x[1], x[1]])
plot(rt.T, ry[:, 0],'x')
show()







