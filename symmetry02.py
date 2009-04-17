from __future__ import division
from numpy import sqrt
from scipy import *
from physicspy.quantum import *


#~ F1 = 1/2
#~ F2 = 3/2
#~ M1 = -1/2
#~ M2 = -3/2
#~ q = M1-M2

#~ # Andy
#~ print "Andy-3j: ", threej(F1,1,F2,-M1,q,M2)
#~ # Cs-text
#~ print "Cs-3j  : ", threej(F2,1,F1,M2,q,-M1)

#~ # Andy
#~ print "Andy-CG: ", sqrt(2*F2+1)*threej(F1,1,F2,-M1,q,M2)
#~ # Cs-text
#~ print "Cs-CG  : ", sqrt(2*F1+1)*threej(F2,1,F1,M2,q,-M1)



#~ F1 = 1/2
#~ M1 = -1/2
#~ F2 = 3/2
#~ M2 = -1/2
#~ q = M2-M1

#~ print (2*F2+1)*threej(F1, 1, F2, M1, q, -M2)**2
#~ print ( 2*F1+1)*threej(F2, 1, F1, M2, -q, -M1)**2


#~ F2 = 5/2
#~ M2 = -5/2


#~ (2*F2+1)*threej(3/2,1,F2,-3/2,(3/2-M2),M2)**2


J2 = 3/2
F2 = 3
M2 = 0

J1 = 1/2
F1  = 4
M1 = arange(-F1,F1+1,1)
s = 0
for mi in M1:
    s += (2*F2+1)*threej(F1,1,F2,-mi,(mi-M2),M2)**2 * (2*F1+1)*(2*J2+1)*sixj(F2, F1, 1, J1, J2, 7/2)**2
print s




F1 = 3/2
M1 = -1/2
F2 = 1/2
M2 = -1/2
#~ print (2*F2+1)*threej(F1,1,F2,-M1,(M1-M2),M2)**2
print threej(F1,1,F2,-M1,(M1-M2),M2)**2
print threej(F2,1,F1,M2,-(M2-M1),-M1)**2


