from __future__ import division
import pylab as pl
from numpy import loadtxt
from threelevel import *


# Excited state lifetimes
# Total should be 1, everything else is scaled with it
G = 1.0
G1 = 0.5 * G
G2 = G - G1
# Rabi frequencies (untis of G)
Om1 = 0.2 * G
Om2 = 0.2 * G
# Laser linewidths (units of G)
gg1 = 0 * G
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

# Comparison
data  = loadtxt("3.txt");
gammaa = 5e6
det = data[:,0]/gammaa

longt = 1000
laser = 1
p1 = sp.zeros((len(det),1))
p2 = sp.zeros((len(det),1))
p3 = sp.zeros((len(det),1))
# Scan the laser detuning
for i in range(0, len(det)):
        if laser == 1 :
            laser1 = (Om1, det[i], gg1)
        else :
            laser2 = (Om2, det[i], gg2)
        params = (G,G1,G2,laser1,laser2)
        temp = sp.dot(expm(Mmat(params) * longt), q0)
        p1[i] = temp[0]
        p2[i] = temp[1]
        p3[i] = temp[2]

pl.figure(1)
pl.plot(det, p3, 'r-', label='P(3)')
pl.plot(det, data[:,1], 'b-', label='Ernie')
pl.xlabel('Detuning (MHz)')
pl.ylabel('Population')
pl.title('Upper state population')

pl.figure(2)
pl.plot(det, p1, 'r-', label='P(1)')
pl.plot(det, data[:,2], 'b-', label='Ernie')
pl.xlabel('Detuning (MHz)')
pl.ylabel('Population')
pl.title('Lower ground state population')

pl.figure(3)
pl.plot(det, p2, 'r-', label='P(2)')
pl.plot(det, data[:,3], 'b-', label='Ernie')
pl.xlabel('Detuning (MHz)')
pl.ylabel('Population')
pl.title('Higher ground state population')

pl.show()
