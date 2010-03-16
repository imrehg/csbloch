#!/usr/bin/env python

from __future__ import division
from numpy import *
from scipy import *
from scipy.linalg import expm
from pylab import *
from scipy.integrate import ode


## Note: Gauss and Square is very different when d >> 0
d12 = 0

Gt = 1
d13 = 0
G31 = 0.9
d23 = 15
G32 = Gt - G31

b = array([[1],[0],[0],[0],[0],[0],[0],[0],[0]])

### Integrate directly

tpulse = 5
npulse = 10

ttotal = tpulse * npulse

om0_1 = 1
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


# Same pulse area (same power (?)) but CW, Square, Gauss
dt1 = 0.05
#rt1, ry1 = doint(b, 0, dt1, ttotal, omconstant)
#rt1, ry1 = doint(b, 0, dt1, ttotal, omsquare)
rt1, ry1 = doint(b, 0, dt1, ttotal, omgauss)
plot(rt1, ry1[:, 0],'b-', linewidth=2, label='g(1)')
plot(rt1, ry1[:, 1],'r-', linewidth=2, label='g(2)')
plot(rt1, ry1[:, 2],'g-', linewidth=2, label='e(3)')
xlim((0, rt1[-1]))

legend(loc = 'best')
show()
