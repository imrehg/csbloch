from pylab import *
from scipy import *
from scipy.linalg import expm, expm2, expm3
from time import time


s = range(3, 100)
t = zeros(size(s))
t2 = zeros(size(s))
t3 = zeros(size(s))
tt = 0
for i in s:
    print "Matrix size %d" %(i)
    m = zeros((i,i))
    half = int(round(i*i/2))+1
    r = random.randint(0,i, (half*2, ))
    n = random.random((half, ))
    for i in range(0, half*2, 2):
        m[r[i], r[i+1]] = n[i/2]
    start = time()
    expm(m)
    t[tt] = time() - start
    start = time()
    expm2(m)
    t2[tt] = time() - start
    start = time()
    expm3(m)
    t3[tt] = time() - start
    tt += 1

print s
print t
plot(s, t, label="Pade")
plot(s, t2, label="Eigenvalue")
plot(s, t3, label="Taylor")
legend(loc = 'best')
show()
