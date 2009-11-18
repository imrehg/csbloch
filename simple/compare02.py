from __future__ import division
import pylab as pl
from numpy import loadtxt
from threelevel import *


# Excited state lifetimes
# Total should be 1, everything else is scaled with it
G = 1.0
G1 = 0.5
G2 = G - G1
# Rabi frequencies (untis of G)
Om1 = 0.2
Om2 = 0.2
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

# Make standard parameters variable
laser1 = (Om1, d1, gg1)
laser2 = (Om2, d2, gg2)
params = (G,G1,G2,laser1,laser2)


(t, p1, p2, p3) = timeseries((0,20,101),params, q0)
pl.plot(t, p1)
pl.plot(t, p2)
pl.plot(t, p3)
pl.show()
