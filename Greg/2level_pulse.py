#!/usr/bin/env python

from __future__ import division
from numpy import *
from scipy import *
from scipy.linalg import expm
from pylab import *
from scipy.integrate import ode


## Note: Gauss and Square is very different when d >> 0
d = 4
G = 1


b = array([[1],[0],[0],[0]])

### Integrate directly

tpulse = 5
npulse = 10

ttotal = tpulse * npulse

om0 = 2
om0t = 0.5

omgc = om0t/10
omgb = om0t/2
omga = om0*om0t / (omgc * sqrt(2*pi))

def omconstant(t):
    ''' Constant Omega '''
    return om0

def omsquare(t):
    ''' Constant Omega '''
    t = mod(t, tpulse)
    if (t < om0t):
        return(om0)
    else:
        return 0

def gauss(x, a, b, c):
    return a * exp(-(x-b)**2/(2*c**2))

def omgauss(t):
    t = mod(t, tpulse)
    return gauss(t, omga, omgb, omgc)



def f(t0, y, omega):
    y1 = zeros(4)
    om = omega(t0)
    obm = array([[0, G, 0, -om], \
                 [0, -G, 0, om], \
                 [0, 0, -1/2*G, d], \
                 [om/2, -om/2, -d, -1/2*G]])
    for i in range(0, 4):
        y1[i] += dot(obm[i,:],y)
    return(y1)


def doint(y0, t0, dt, tend, omega):
    r = ode(f).set_integrator('zvode', with_jacobian=False)
    r.set_initial_value(y0, t0).set_f_params(omega)
    rt = array([t0])
    ry = array(y0.T)
    while r.successful() and r.t < tend:
        r.integrate(r.t+dt)
        rt = append(rt, r.t)
        ry = append(ry, r.y.T, 0)
    return (rt.T, ry)



#dt1 = 0.01
#rt1, ry1 = doint(b, 0, dt1, 5, omconstant)
#plot(rt1, ry1[:, 0],'bx', label='Ground state, dt = %.2f'%(dt1))
#plot(rt1, ry1[:, 1],'ro', label='Upper state, dt = %.2f'%(dt1))


dt1 = 0.01
rt1, ry1 = doint(b, 0, dt1, ttotal, omgauss)
plot(rt1, ry1[:, 0],'b-', label='Ground state, dt = %.2f'%(dt1))
plot(rt1, ry1[:, 1],'b--', label='Upper state, dt = %.2f'%(dt1))

dt2 = 0.01
rt2, ry2 = doint(b, 0, dt2, ttotal, omsquare)
plot(rt2, ry2[:, 0],'r-', label='Ground state, dt = %.2f'%(dt2))
plot(rt2, ry2[:, 1],'r--', label='Upper state, dt = %.2f'%(dt2))
#legend(loc = 'best')


#t = linspace(0, 5, 1000)
#y = omgauss(t)
#plot(t, y)
show()
