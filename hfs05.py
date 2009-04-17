#!/usr/bin/env python
from __future__ import division
from numpy import *
from pylab import *
from scipy import *
from scipy import sparse
from scipy.linalg import expm
from physicspy.quantum import *

def Ehf(F,I,J,A):
	return A/2*(F*(F+1) - I*(I+1) - J*(J+1))

def gf(F,J,I):
    return (F*(F+1) + J*(J+1) - I*(I+1))/(2*F*(F+1))

def i2n(i,j,n):
    return i*n + j
    

def n2i(ni, n):
    return array([floor(ni/n), (ni%n)])




I = 7/2
gI = -0.00039885395
gS = 2.0023193043622
gL = 0.99999587
uB = 1.399624e6
#~ J = array([1/2, 1/2, 3/2])
J = array([1/2, 1/2])
#~ J = array([1/2])
Ahfs = array([2.2981572425e9,  291.9201e6, 50.28827e6])
Bhfs = array([0, 0, -0.4934e6])
Chfs = array([0, 0, 0.560e3])
gJ = array([2.00254032, 0.665900, 1.334])
baseE = array([0, 0, 0])
names  = array(['S1/2', 'P1/2', 'P3/2'])
B = 1

lev = array([])
print "J....F....mf....dEhfs(MHz)....dE(MHz)/G"
for i in range (0,J.size):
    #~ print "--------------"
    #~ print names[i]," hyperfine states"
    F = arange(abs(I-J[i]),I+J[i]+1,1)
    #~ print "F   = ",F
    #~ energy = Ehf(F,I,J[i],Ahfs[i])
    #~ print "Ehf = ",energy
    for fi in range(0,F.size):
        energy = energyhfsalkali(F[fi],J[i],I,Ahfs[i],Bhfs[i],Chfs[i])
        gf = landeg(gJ[i],gI,F[fi],I,J[i])
        EB = uB*gf
        mf = arange(-F[fi],F[fi]+1,1)
        for mi in mf:
            lev = append(lev,[i,J[i], F[fi], mi, energy+baseE[i]+EB*mi*B])


lev = reshape(lev,[-1,5])
#~ print lev[:,3:]


#~ t = 
#~ A21 = 1/t/2 
n = lev.shape[0]
print "Number of sublevels: ",n
A = sparse.lil_matrix((n**2,n**2),dtype=complex)

ng = 16
ne = 16

q = 1
Om = 1
detu = +4514.7351
lev[:,4] /= 1e6

#~ for j in range(0,n):
    #~ for k in range(0,n):
        #~ delta = lev[k,4] - lev[j,4]
        #~ A[i2n(j,k,n),i2n(j,k,n)] = -1j*delta
        #~ for l in range(0,n):
            #~ if (int(lev[j,0]) != (lev[l,0])):
                #~ if 
            
            
            
            
            
#~ print A


### Sorting hyperfine energy sub-levels within Zeeman-levels
x = lev[0:16,:]
xs = argsort(x[:,4])
lev[0:16,:] = x[xs,:]
x = lev[16:32,:]
xs = argsort(x[:,4])
lev[16:32,:] = x[xs,:]

#~ print lev[:,3:]

def ddetu(k,j):
    delta = detu - (lev[k,4] - lev[j,4])
    if abs(delta) > 300 :
        return 0
    else:
        return 1



for j in range(0,ng):
    for k in range(ng,n):
        delta = detu - (lev[k,4] - lev[j,4])
        A[i2n(j,k,n),i2n(j,k,n)] = -1j*delta
        A[i2n(k,j,n),i2n(k,j,n)] = 1j*delta
        
for j in range(0,n):
    for k in range(0,n):
        for l in range(0,n):
            if (int(lev[j,0]) != (lev[l,0])):
                if ((lev[max([j,l]),3] - lev[min([j,l]),3]) == q) and (abs(lev[j,2] - lev[l,2])<2):
                    A[i2n(j,k,n),i2n(l,k,n)]  += 1j/2*Om * ddetu(j,l) 
            if (int(lev[l,0]) != (lev[k,0])):
                if ((lev[max([l,k]),3] - lev[min([l,k]),3]) == q) and (abs(lev[l,2] - lev[k,2])<2):
                    A[i2n(j,k,n),i2n(j,l,n)]  -= 1j/2*Om * ddetu(k,l)


A[0,:] = zeros(n**2)
for j in range(0,n):
    A[0,i2n(j,j,n)] = 1
    
print A

obm = A.todense()
b = zeros(n**2);
j0 = 0
b[i2n(j0,j0,n)] = 1

t = linspace(0,1e-6,3)
p1 = []
p2 = []
for k in range(0,t.size):
    print k
    obmt = expm(obm*t[k])
    p1.append(dot(obmt,b)[0])
p1 = array(p1)
plot(t,p1)
#~ plot(t,p2,[t.min(), t.max()],[x[1], x[1]])
show()
