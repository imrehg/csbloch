#!/usr/bin/env python
"""
Reproducing results from:
Soares and Araujo, JOSA:B 43 (2010) 085003
Autler-Townes doublet and electromagnetically induced transparency resonance
probed by an ultrashort pulse train

For testing purposes, to see whether our algorithm works well.
The level structure slightly different (i.e. ordered in increasing energy)

    ---- (3)
    /  \
   /   ----(2)
----(1) 

Compared to the paper: (abc) -> (312)

"""

from __future__ import division
from pysparse.spmatrix import *
import numpy as np
import physicspy.quantum as quantum

def blochlambda(G, G12=0, gl=[0, 0]):
    """
    genmat(G)
    Generate the parts of the Bloch-equations to use in the
    actual calculation:
    use atomic levels, but ignore: laser linewidth, ground state decoherence

    Input:
    G = [G31, G32] the upper level decay rate to the two ground states
    G12 = decoherence between the ground states (default = 0)
    gl = laser linewidth for laser 13 and 23 (default = [0 0])


    Output:
    bloch_atomic, laser1_pow, laser1_det, laser2_pow, laser2_det, laser12_det
    """
    G31 = G[0]
    G32 = G[1]
    Gtot = G31 + G32

    bloch_atomic = ll_mat(9, 9, 10)
    bloch_atomic[0, 2] = G31
    bloch_atomic[1, 2] = G32
    bloch_atomic[2, 2] = -Gtot
    cdecay = [4, 5, 7, 8]
    for i in cdecay:
        bloch_atomic[i, i] = -Gtot / 2 
    # Ground state decoherence
    bloch_atomic[3, 3] += -G12 / 2
    bloch_atomic[6, 6] += -G12 / 2
    # Laser linewidths
    bloch_atomic[3, 3] += -(gl[0] + gl[1]) / 2
    bloch_atomic[4, 4] += -gl[0] / 2
    bloch_atomic[5, 5] += -gl[1] / 2
    bloch_atomic[6, 6] += -(gl[0] + gl[1]) / 2
    bloch_atomic[7, 7] += -gl[0] / 2
    bloch_atomic[8, 8] += -gl[1] / 2


    laser1_pow = ll_mat(9, 9, 10)
    laser1_pow[0, 7] = 1
    laser1_pow[2, 7] = -1
    laser1_pow[3, 8] = 1/2
    laser1_pow[5, 6] = -1/2
    laser1_pow[6, 5] = 1/2
    laser1_pow[7, 0] = -1/2
    laser1_pow[7, 2] = 1/2
    laser1_pow[8, 3] = -1/2
    laser1_det = ll_mat(9, 9, 10)
    laser1_det[4, 7] = 1
    laser1_det[7, 4] = -1

    laser2_pow = ll_mat(9, 9, 10)
    laser2_pow[1, 8] = 1
    laser2_pow[2, 8] = -1
    laser2_pow[3, 7] = 1/2
    laser2_pow[4, 6] = 1/2
    laser2_pow[6, 4] = -1/2
    laser2_pow[7, 3] = -1/2
    laser2_pow[8, 1] = -1/2
    laser2_pow[8, 2] = 1/2
    laser2_det = ll_mat(9, 9, 10)
    laser2_det[5, 8] = 1
    laser2_det[8, 5] = -1

    laser12_det = ll_mat(9, 9, 10)
    laser12_det[3, 6] = 1
    laser12_det[6, 3] = -1

    return bloch_atomic, laser1_pow, laser1_det, laser2_pow, laser2_det, laser12_det


def blochd1():
    """
    Cesium D1 line @ B=0
    Simulation reduces to 4 levels F={3,4} -> F'={3,4}
    """

    G = [0, 0, 1, 1]
    Jlower = 1/2
    Jupper = 1/2
    I = 7/2
    Flower = [3, 4]
    Fupper = [3, 4]
    nlev = len(Flower) + len(Fupper)
    # Decay matrix with branching
    A = np.zeros((nlev, nlev))
    off = len(Flower)
    # branching ratio
    for i in xrange(len(Fupper)):
        for f in xrange(len(Flower)):
            prob = (2*Jupper+1)*(2*Flower[f]+1)*quantum.sixj(Jupper, Jlower, 1, Flower[f], Fupper[i], I)**2
            A[f, i+off] = prob

    M = np.zeros((nlev, nlev))
    off = len(Flower)
    # branching ratio
    for i in xrange(len(Flower)):
        for f in xrange(len(Fupper)):
            prob = np.sqrt((2*Jlower+1)*(2*Fupper[f]+1)*quantum.sixj(Jlower, Jupper, 1, Fupper[f], Flower[i], I)**2)
            M[i, f+off] = prob
            M[f+off, i] = prob

    # ### Temporary
    # for i in xrange(4):
    #     M[i,3] = 0
    #     M[3,i] = 0
    # for i in xrange(4):
    #     for j in xrange(4):
    #         M[i, j] = np.sign(M[i, j])*0.5
    # ###

    # The diagonal states
    # Helpers
    dg = []
    xlook = {}
    ylook = {}
    ymulti = {}
    n = 0
    for i in xrange(nlev-1):
        for j in xrange(i+1, nlev):
            dg += [[i, j]]
            xlook["%d:%d"%(i,j)] = n
            xlook["%d:%d"%(j,i)] = n
            ylook["%d:%d"%(i,j)] = n
            ylook["%d:%d"%(j,i)] = n
            n += 1
    ldg = len(dg)

    bloch_atomic = ll_mat(16, 16, 20)
    for i in xrange(nlev):
        # ### Temporary:
        # if i == 3:
        #     continue
        # ###
        bloch_atomic[i, i] += -G[i]
    for i in xrange(nlev):
        for j in xrange(nlev):
            # ### Temporary:
            # if (i == 3) | (j == 3):
            #     continue
            # ###
            bloch_atomic[j, i] += A[j, i]

    for i in xrange(ldg):
        # ### Temporary:
        # if i in [2, 3, 5]:
        #     continue
        # ###
        bloch_atomic[i+nlev, i+nlev] += -(G[dg[i][0]] + G[dg[i][1]])/2
        bloch_atomic[i+nlev+ldg, i+nlev+ldg] += -(G[dg[i][0]] + G[dg[i][1]])/2

    laser_pow = ll_mat(16, 16, 20)
    for j in xrange(nlev):
        for l in xrange(nlev):
            if l != j:
                ypos = nlev + ldg + ylook["%d:%d"%(l, j)]
                mul = 1 if l>j else -1
                laser_pow[j, ypos] += mul*2*M[l, j]

    for j in xrange(nlev-1):
        for k in xrange(j+1, nlev):
            xjk = nlev + xlook["%d:%d"%(j, k)]
            for l in xrange(nlev):
                if l !=k:
                    ylk = nlev + ldg + ylook["%d:%d" %(l, k)]
                    scale = 1 if l < k else -1
                    laser_pow[xjk, ylk] += -M[j, l]*scale
                if j !=l:
                    yjl = nlev + ldg + ylook["%d:%d" %(j, l)]
                    scale = 1 if j < l else -1
                    laser_pow[xjk, yjl] += M[l, k]*scale

    for j in xrange(nlev-1):
        for k in xrange(j+1, nlev):
            yjk = nlev + ldg + ylook["%d:%d"%(j, k)]
            for l in xrange(nlev):
                if l !=k:
                    xlk = nlev + xlook["%d:%d" %(l, k)]
                    laser_pow[yjk, xlk] += M[j, l]
                else:
                    laser_pow[yjk, l] += M[j, l]

                if j !=l:
                    xjl = nlev + xlook["%d:%d" %(j, l)]
                    laser_pow[yjk, xjl] += -M[l, k]
                else:
                    laser_pow[yjk, l] += -M[l, k]

    return (bloch_atomic, laser_pow,)

def laser_det(d):
    nlev = 4
    ldg = 3*2
    laser_det = ll_mat(16, 16, 5)
    ## 01, 02, 03, 12, 13, 23
    d01 = -2009.3
    d23 = -255.23
    # detunings = np.array([d01, d, d+d23, d-d01, d-d01+d23, d23])*2*np.pi
    detunings = np.array([d01, d, d+d23, d-d01, d-d01+d23, d23])
    for i in xrange(ldg):
        laser_det[i+nlev, i+nlev+ldg] = detunings[i]
        laser_det[i+nlev+ldg, i+nlev] = -detunings[i]
    return laser_det


def blochd2():
    """
    Cesium D1 line @ B=0
    Simulation reduces to 4 levels F={3,4} -> F'={3,4}
    """

    G = [0, 0, 1, 1, 1, 1]
    Jlower = 1/2
    Jupper = 3/2
    I = 7/2
    Flower = [3, 4]
    Fupper = [2, 3, 4, 5]
    nlev = len(Flower) + len(Fupper)
    # Decay matrix with branching
    A = np.zeros((nlev, nlev))
    off = len(Flower)
    # branching ratio
    for i in xrange(len(Fupper)):
        for f in xrange(len(Flower)):
            prob = (2*Jupper+1)*(2*Flower[f]+1)*quantum.sixj(Jupper, Jlower, 1, Flower[f], Fupper[i], I)**2
            A[f, i+off] = prob

    M = np.zeros((nlev, nlev))
    off = len(Flower)
    # branching ratio
    for i in xrange(len(Flower)):
        for f in xrange(len(Fupper)):
            prob = np.sqrt((2*Jlower+1)*(2*Fupper[f]+1)*quantum.sixj(Jlower, Jupper, 1, Fupper[f], Flower[i], I)**2)
            M[i, f+off] = prob
            M[f+off, i] = prob

    # ### Temporary
    # for i in xrange(nlev):
    #     for j in xrange(nlev):
    #         if (i <> 3) & (j <> 3):
    #             M[i,j] = 0
    #             M[j, i] = 0
    # print M
    # for i in xrange(4):
    #     for j in xrange(4):
    #         M[i, j] = np.sign(M[i, j])*0.5
    # ###

    # The diagonal states
    # Helpers
    dg = []
    xlook = {}
    ylook = {}
    ymulti = {}
    n = 0
    for i in xrange(nlev-1):
        for j in xrange(i+1, nlev):
            dg += [[i, j]]
            xlook["%d:%d"%(i,j)] = n
            xlook["%d:%d"%(j,i)] = n
            ylook["%d:%d"%(i,j)] = n
            ylook["%d:%d"%(j,i)] = n
            n += 1
    ldg = len(dg)

    bloch_atomic = ll_mat(36, 36, 20)
    for i in xrange(nlev):
        # ### Temporary:
        # if i == 3:
        #     continue
        # ###
        bloch_atomic[i, i] += -G[i]
    for i in xrange(nlev):
        for j in xrange(nlev):
            # ### Temporary:
            # if (i == 3) | (j == 3):
            #     continue
            # ###
            bloch_atomic[j, i] += A[j, i]

    for i in xrange(ldg):
        # ### Temporary:
        # if i in [2, 3, 5]:
        #     continue
        # ###
        bloch_atomic[i+nlev, i+nlev] += -(G[dg[i][0]] + G[dg[i][1]])/2
        bloch_atomic[i+nlev+ldg, i+nlev+ldg] += -(G[dg[i][0]] + G[dg[i][1]])/2

    laser_pow = ll_mat(36, 36, 20)
    for j in xrange(nlev):
        for l in xrange(nlev):
            if l != j:
                ypos = nlev + ldg + ylook["%d:%d"%(l, j)]
                mul = 1 if l>j else -1
                laser_pow[j, ypos] += mul*2*M[l, j]

    for j in xrange(nlev-1):
        for k in xrange(j+1, nlev):
            xjk = nlev + xlook["%d:%d"%(j, k)]
            for l in xrange(nlev):
                if l !=k:
                    ylk = nlev + ldg + ylook["%d:%d" %(l, k)]
                    scale = 1 if l < k else -1
                    laser_pow[xjk, ylk] += -M[j, l]*scale
                if j !=l:
                    yjl = nlev + ldg + ylook["%d:%d" %(j, l)]
                    scale = 1 if j < l else -1
                    laser_pow[xjk, yjl] += M[l, k]*scale

    for j in xrange(nlev-1):
        for k in xrange(j+1, nlev):
            yjk = nlev + ldg + ylook["%d:%d"%(j, k)]
            for l in xrange(nlev):
                if l !=k:
                    xlk = nlev + xlook["%d:%d" %(l, k)]
                    laser_pow[yjk, xlk] += M[j, l]
                else:
                    laser_pow[yjk, l] += M[j, l]

                if j !=l:
                    xjl = nlev + xlook["%d:%d" %(j, l)]
                    laser_pow[yjk, xjl] += -M[l, k]
                else:
                    laser_pow[yjk, l] += -M[l, k]

    return (bloch_atomic, laser_pow,)

def laser_detd2(d):
    dlow = [0, 1756.33]
    dhigh = [0, 28.89, 67.35, 115.32]
    total = [-dl for dl in dlow]+[d-dh for dh in dhigh]
    nlev = len(total);
    ncross = int(nlev*(nlev-1)/2)
    detunings = []
    for i in xrange(nlev-1):
        for j in xrange(i+1, nlev):
            detunings += [total[j] - total[i]]
    detu_matrix = ll_mat(nlev*nlev, nlev*nlev, nlev*2)
    for i in xrange(ncross):
        detu_matrix[i+nlev, i+nlev+ncross] = detunings[i]
        detu_matrix[i+nlev+ncross, i+nlev] = -detunings[i]

    return detu_matrix
