from __future__ import division
import pylab as pl
from numpy import loadtxt
from threelevelt import *
import scipy as sp
from scipy.integrate import odeint

# Excited state lifetimes
# Total should be 1, everything else is scaled with it
G = 1.0
G1 = 0.3
G2 = G - G1
# Rabi frequencies (untis of G)
Om1 = 10000
Om2 = Om1
# Laser linewidths (units of G)
gg1 = 0
gg2 = gg1
# Detuning (Units of G)
d1 = 0
d2 = -1
# Starting density matrix elements
# Matrix elements: x11, x22, x33, x12, x13, x23, y12, y13, y23
# x11, x22, x33 (meaning: q0[0:3]) are the populations
q0 = sp.zeros((9,1))
# Start in level 1
q0[0] = 1

dg = 1838
d1 = -dg/2

# Make standard parameters variable
# laser1 = (Om1, Om2, d1, gg1)
# params = (G,G1,G2,dg,laser1)

# M = MmatSpecial(params)

# (t, p1, p2, p3) = timeseries((0,20,101),params, q0)
# pl.plot(t, p1)
# pl.plot(t, p2)
# pl.plot(t, p3)
# pl.show()

laser1 = (Om1, Om2, d1, gg1)
params = (G,G1,G2,dg,laser1)
Mzero, Mom = MmatSpecial(params)

ev = lambda t, tg: sp.exp(-(t/tg)**2)

def dypulse(y, t, params):
    Mzero, Mom, Om, tg = params
    M = Mzero + ev(t, tg)*Om*Mom
    return sp.dot(M, y.T)

tp = 1/200000
tg = tp/sp.sqrt(2*sp.log(2))
tl = tp*4
tsteps = sp.linspace(-tl, tl, 1000)
trep = 1/18.38
t = sp.linspace(0, trep, 11)

y0 = sp.zeros((9))
y0[0] = 1/2
y0[1] = 1/2
y0[2] = 0
nsteps = 600*2

tt = sp.zeros((0))
yo0 = sp.zeros((0))
yo1 = sp.zeros((0))
yo2 = sp.zeros((0))
for n in xrange(nsteps):
    y = odeint(dypulse, y0, tsteps, args=((Mzero,Mom,Om1,tg),))
    y0 = y[-1]
    y = odeint(dypulse, y0, t, args=((Mzero,Mom,0,10),))
    tt = sp.append(tt, t+n*trep)
    yo0 = sp.append(yo0, y.T[0])
    yo1 = sp.append(yo1, y.T[1])
    yo2 = sp.append(yo2, y.T[2])
    y0 = y[-1]

pl.plot(tt, yo0)
pl.plot(tt, yo1)
pl.plot(tt, yo2)
pl.show()

# pl.plot(tsteps, ev(tsteps, tg))
# pl.show()





# laser1 = (0, 0, d1, gg1)
# params = (G,G1,G2,dg,laser1)
# Moff = MmatSpecial(params)

# t = sp.linspace(0,5,8001)

# y = odeint(dy, y0, t, args=(Mon,))
# yt = y.T
# # y2 = []
# # for yc in y:
# #     y2.append(yc[2])
# pl.plot(t, yt[0])
# pl.plot(t, yt[1])
# pl.plot(t, yt[2])

# y = odeint(dy, y0, t, args=(Moff,))
# yt = y.T
# # y2 = []
# # for yc in y:
# #     y2.append(yc[2])
# pl.plot(t, yt[0])
# pl.plot(t, yt[1])
# pl.plot(t, yt[2])


# pl.show()
