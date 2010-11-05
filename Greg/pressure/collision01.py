from __future__ import division
import numpy as np
import pylab as pl
from scipy.linalg import expm, pinv, lstsq

"""
Trying to recreate and explore the theory in
Sung and Bergman, PRA 39 (1989) 6298
Theory of optical coherent transients including collisional
effects: application to an extended-pulse photon echo
"""

def blochmatrix(G, Om, det, Gcoll=0, Gphase=0):
    bloch = np.zeros((4,4))
    # Parameter order: x11, x22, x12, y12
    bloch[0, 0] = -Gcoll
    bloch[0, 1] = G
    bloch[0, 3] = 2 * Om

    bloch[1, 1] = -(G + Gcoll)
    bloch[1, 3] = -2 * Om

    bloch[2, 2] = -(G/2 + Gcoll + Gphase)
    bloch[2, 3] = det

    bloch[3, 0] = -Om
    bloch[3, 1] = Om
    bloch[3, 2] = -det
    bloch[3, 3] = -(G/2 + Gcoll + Gphase)
    return bloch


def getfinalpop(*args, **kargs):
    bm = blochmatrix(*args, **kargs)    
    zz = np.zeros((5,1))
    zz[-1] = 1
    mm = np.zeros((5, 4))
    for i in xrange(4):
        for j in xrange(4):
            mm[i, j] = bm[i,j]
    for j in xrange(2):
        mm[4, j] = 1
    final = lstsq(mm, zz)
    return final[0][0:2]

def detuv(v):
    bmv = np.zeros((4,4))
    bmv[2, 3] += v
    bmv[3, 2] += -v
    return bmv
    
def wv(v1, v2, Gcoll, sigma):
    return Gcoll/sigma/np.sqrt(np.pi)*np.exp(-(v1-v2)**2/sigma**2)

def calcchange(v, bmv, rhov, Gcoll, sigma, dt):
    
    dv = v[1] - v[0]
    for i in xrange(len(v)):
        ww = 0
        for j in xrange(len(v)):
            wj = wv(v[i], v[j], Gcoll, sigma) 
            for k in xrange(4):
                bmv[i][k,k] += wj*rhov[j][k]*dv
    print bmv[1]

def collisiontime(G, Om, det, Gcoll, Gphase, v, vp, sigma, t):
    bm = blochmatrix(G, Om, det, Gcoll, Gphase)
    bmv = [bm+detuv(vi) for vi in v]
    dv = v[1]-v[0]

    rhov = []
    for i, vi in enumerate(v):
        start = np.zeros((4,1))
        start[0, 0] = vp[i]
        rhov += [start]

    calcchange(v, bmv, rhov, Gcoll, sigma, 0)
    # for i in xrange(len(t)):

    # print wv(v[0], v[1], Gcoll, sigma)*dv
    



if __name__ == "__main__":

    # ## Phase changing collisions
    # detlist = np.linspace(-35, 35, 2001)
    # G = 1
    # Om = 2
    # pop3 = []
    # for det in detlist:
    #     pop3 += [getfinalpop(G, Om, det)[1]]
    # pl.plot(detlist, pop3)
    # print sum(pop3)*(detlist[1]-detlist[0])

    # G = 1
    # Om = 1
    # pop3 = []
    # for det in detlist:
    #     pop3 += [getfinalpop(G, Om, det, Gcoll=0, Gphase=1)[1]]
    # pl.plot(detlist, pop3)
    # pl.show()

    G = 1
    Om = 1
    det = 0
    vv = np.linspace(-5, 5, 3)
    vp = np.zeros(len(vv))
    vp[(len(vv)-1)/2] = 1
    collisiontime(G, Om, det, 1, 0, vv, vp, 1, 1)






    # bm = blochmatrix(G, Om, det)
    # start = np.zeros((4,1))
    # start[0, 0] = 1

    
    # zz = np.zeros((5,1))
    # zz[-1] = 1
    # mm = np.zeros((5, 4))
    # for i in xrange(4):
    #     for j in xrange(4):
    #         mm[i, j] = bm[i,j]
    # for j in xrange(2):
    #     mm[4, j] = 1
    # final = lstsq(mm, zz)
    # print final[0][0]
    # tlist = np.linspace(0, 10, 1001)
    # res = [np.dot(expm(bm*t), start) for t in tlist]
    # pop1 = []
    # pop2 = []
    # for r in res:
    #     pop1 += [r[0]]
    #     pop2 += [r[1]]
    # pl.plot(tlist, pop1)
    # pl.plot(tlist, pop2)
    # pl.plot([tlist[0], tlist[-1]], [final[0][0], final[0][0]])
    # pl.plot([tlist[0], tlist[-1]], [final[0][1], final[0][1]])
    # pl.show()



    
