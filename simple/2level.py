import scipy as sp
from scipy.linalg import *
import pylab as pl

# Excited state lifetime
# Should be 1, everything else is scaled with it
G = 1
# Rabi frequency
Om1 = 2.2
# Laser linewidth
gg = 0.1
# Detuning
d = -0.5

M = sp.array([[0,1j*(-Om1/2),1j*(Om1/2),G],
              [1j*(-Om1/2),-(G+gg)/2+1j*d,0,1j*(Om1/2)],
              [1j*(Om1/2),0,-(G+gg)/2-1j*d,1j*(-Om1/2)],
              [0,1j*(Om1/2),1j*(-Om1/2),-G]])

q0 = sp.array([[0],[0],[0],[1]])
t = sp.linspace(0,10,101)
p1 = sp.zeros((len(t),1))
p2 = sp.zeros((len(t),1))

for i in range(0, len(t)):
    temp = sp.dot(expm(M * t[i]), q0)
    p1[i] = temp[0]
    p2[i] = temp[3]

pl.plot(t,p1,'b--',label='Ground state')
pl.plot(t,p2,'r-',label='Excited state')
pl.xlabel('time')
pl.ylabel('Population')
pl.ylim([0,1])
pl.legend(loc="lower right")
pl.show()
