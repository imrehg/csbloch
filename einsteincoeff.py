from numpy import *
from physicspy import *

c = SPEEDOFLIGHT()
hbar = HBAR()
e0 = 8.854187817e-12

l = 894e-9
D12 = 2.7e-29
J = 0.5

l2 = 852e-9
D122 = 3.8e-29
J2 = 1.5


def Acoeff(l,D12,J):
    w = c/l*2*pi
    return w**3*2/(2*J+1)*D12**2/(3*pi*e0*hbar*c**3)


A = Acoeff(l,D12,J)
i01 = 8 * pi ** 3 * hbar * (A / (2 * pi)) * c / (3 * l ** 3)

A2 = Acoeff(l2,D122,J2)
i02 = 8 * pi ** 3 * hbar * (A2 / (2 * pi)) * c / (3 * l2 ** 3)
print i01, i02




l = 397e-9
t = 7.098e-9
p = 0.95

#~ l = 866e-9
#~ t = 7.098e-9
#~ p = 0.05

G = 1/t/(2*pi)
I_0 = 4*pi**2*hbar*c/(3*t*l**3)
r = 6*pi*c**3*(1/t*p)/(hbar*(c/l*2*pi)**3)*I_0/c

print I_0, sqrt(r)/(2*pi*1e6)