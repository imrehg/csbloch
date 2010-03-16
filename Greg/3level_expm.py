#!/usr/bin/env python

from __future__ import division
from numpy import *
from scipy import *
from scipy.linalg import expm
from pylab import *
from scipy.integrate import ode
from time import time

## Note: Gauss and Square is very different when d >> 0
d12 = 0

Gt = 1
d13 = 0
G31 = 0.9
d23 = 0
G32 = Gt - G31

b = array([[1],[0],[0],[0],[0],[0],[0],[0],[0]])

### Integrate directly

tpulse = 5
npulse = 5

ttotal = tpulse * npulse

om0_1 = 2
om0_2 = 1
om0t = 1

omgc = om0t/20
omgb = om0t/2
omga = om0_1*om0t / (omgc * sqrt(2*pi))

om0_1c = om0_1 * om0t / tpulse
om0_2c = om0_2 * om0t / tpulse

def omconstant(t):
    ''' Constant Omega '''
    return (om0_1c, om0_2c)

def omsquare(t):
    ''' Constant Omega '''
    t = mod(t, tpulse)
    if (t < om0t):
        return(om0_1, om0_1)
    else:
        return(0, 0)

def gauss(x, a, b, c):
    return a * exp(-(x-b)**2/(2*c**2))

def omgauss(t):
    t = mod(t, tpulse)
    om = gauss(t, omga, omgb, omgc)
    return (om, om)

def f(t0, y, omega):
    y1 = zeros(9)
    om13, om23 = omega(t0)
    Gtot = G31 + G32
    # meaning of d12?
    d12 = 0
    # Components of the state vector:
    # x11, x22, x33, x12, x13, x23, y12, y13, y23
    # Non-zer0: om31, om32
    obm = array([[0, 0, G31, 0, 0, 0, 0, om13, 0], \
                 [0, 0, G32, 0, 0, 0, 0, 0, om23], \
                 [0, 0, -Gtot, 0, 0, 0, 0, -om13, -om23], \
                 [0, 0, 0, 0, 0, 0, d12, om23/2, om13/2], \
                 [0, 0, 0, 0, -Gtot/2, 0, om23/2, d13, 0], \
                 [0, 0, 0, 0, 0, -Gtot/2, -om13/2, 0, d23], \
                 [0, 0, 0, -d12, -om23/2, om13/2, -Gtot/2, 0, 0], \
                 [-om13/2, 0, om13/2, -om23/2, -d13, 0, 0, -Gtot/2, 0], \
                 [0, -om23/2, om23/2, -om13/2, 0, -d23, 0, 0, -Gtot/2]])
    for i in range(0, 9):
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


def fmat(t0, dt, omega):
    y1 = zeros(9)
    om13, om23 = omega(t0+dt/2)
    Gtot = G31 + G32
    # meaning of d12?
    d12 = 0
    # Components of the state vector:
    # x11, x22, x33, x12, x13, x23, y12, y13, y23
    # Non-zer0: om31, om32
    obm = array([[0, 0, G31, 0, 0, 0, 0, om13, 0], \
                 [0, 0, G32, 0, 0, 0, 0, 0, om23], \
                 [0, 0, -Gtot, 0, 0, 0, 0, -om13, -om23], \
                 [0, 0, 0, 0, 0, 0, d12, om23/2, om13/2], \
                 [0, 0, 0, 0, -Gtot/2, 0, om23/2, d13, 0], \
                 [0, 0, 0, 0, 0, -Gtot/2, -om13/2, 0, d23], \
                 [0, 0, 0, -d12, -om23/2, om13/2, -Gtot/2, 0, 0], \
                 [-om13/2, 0, om13/2, -om23/2, -d13, 0, 0, -Gtot/2, 0], \
                 [0, -om23/2, om23/2, -om13/2, 0, -d23, 0, 0, -Gtot/2]])
    obmt = expm(obm * dt)
#    for i in range(0, 9):
#        y1[i] += dot(obmt[i,:],y)
#    return(y1, obmt)
    return(obmt)


def solveexp(y0, tpulse, npulse, omega):
    obm0 = identity(9)
    pt = linspace(0, om0t, 41)
    dt = pt[1] - pt[0]
    for t in pt[0:-1]:
        obm0 = dot(fmat(t,dt, omega), obm0)
    obm0 = dot(fmat(om0t, tpulse-om0t, omega), obm0)
    y1 = y0.T
    x = zeros((1, 9))
    for i in range(0, 9):
        x[0,i] = dot(obm0[i,:],y1.T)
    y1 = append(y1, x, 0)
    for n in range(1, npulse):
        x = zeros((1, 9))
        for i in range(0, 9):
            x[0,i] = dot(obm0[i,:],y1[n,:].T)
        y1 = append(y1, x, 0)
    t1 = arange(0, npulse+1) * tpulse
    return(t1.T, y1)


# Same pulse area (same power (?)) but CW, Square, Gauss
dt1 = 0.05
start = time()
#rt1, ry1 = doint(b, 0, dt1, ttotal, omconstant)
#rt1, ry1 = doint(b, 0, dt1, ttotal, omsquare)
rt1, ry1 = doint(b, 0, dt1, ttotal, omgauss)
print "Full-numeric: %f s" %(time() - start)


start = time()
#te, ye = solveexp(b, tpulse, npulse, omsquare)
te, ye = solveexp(b, tpulse, npulse, omgauss)
print "Matrix exponential: %f s" %(time() - start)



plot(rt1, ry1[:, 0],'b-', linewidth=2, label='g(1)')
plot(rt1, ry1[:, 1],'r-', linewidth=2, label='g(2)')
plot(rt1, ry1[:, 2],'g-', linewidth=2, label='e(3)')
xlim((0, rt1[-1]))
for i in range(0, 3):
    plot(te, ye[:, i], 's')

legend(loc = 'best')
show()
