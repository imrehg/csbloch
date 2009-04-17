from numpy import sqrt
from physicspy.quantum import *


F1 = 1/2
F2 = 3/2
M1 = -1/2
M2 = -3/2
q = M1-M2

# Andy
print "Andy-3j: ", threej(F1,1,F2,-M1,q,M2)
# Cs-text
print "Cs-3j  : ", threej(F2,1,F1,M2,q,-M1)

# Andy
print "Andy-CG: ", sqrt(2*F2+1)*threej(F1,1,F2,-M1,q,M2)
# Cs-text
print "Cs-CG  : ", sqrt(2*F1+1)*threej(F2,1,F1,M2,q,-M1)
