from numpy import zeros, float64
from pylab import *
from scipy import random
from scipy.linalg import expm, expm2, expm3
from expf import expf
from time import time


s = range(3, 100)
t = zeros(size(s))
t2 = zeros(size(s))
t3 = zeros(size(s))
t4 = zeros(size(s))
tt = 0
for i in s:
    print "Matrix size %d" %(i)
    m = zeros((i,i), float64)
    half = int(round(i*i/2))+1
    r = random.randint(0,i, (half*2, ))
    n = random.random((half, ))
    for z in range(0, half*2, 2):
        m[r[z], r[z+1]] = n[z/2]
    start = time()
    expm(m)
    t[tt] = time() - start

    start = time()
    expm2(m)
    t2[tt] = time() - start

    start = time()
    expm3(m)
    t3[tt] = time() - start

    start = time()
    expf(i, 1,  m, i)
    t4[tt] = time() - start
    tt += 1

print s
print t
plot(s, t, label="Pade")
plot(s, t2, label="Eigenvalue")
plot(s, t3, label="Taylor")
plot(s, t4, label="Fortran")
legend(loc = 'best')
show()
