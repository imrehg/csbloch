#!/usr/bin/env python
"""
Reproducing results from:
Matthew McDonnell's thesis, 2003
http://www.physics.ox.ac.uk/Users/iontrap/publications.html

    ---- (3)
    /  \
   /   ----(2)
----(1) 


"""

from __future__ import division
from numpy import *
from scipy import *
from scipy.linalg import expm
import pylab as pl
from scipy.integrate import ode, quad
from time import time, clock
from pysparse.spmatrix import *
from expv import spexpv
from blochmatrix import *

def timedomain(args):
    """
    Input:
    mnames, t, v0, det13, det23, omega13, omega23

    """
    mnames, t, v0, det13, det23, omega13, omega23 = args
    for name in mnames:
        globals()[name] = ll_mat_from_mtx(name)
    # detunings
    d13real = det13
    d23real = det23
    d12real = d13real - d23real

    # Between pulses: CW laser is on
    bloch = bloch_atomic.copy()
    bloch.shift(omega13, laser1_pow)
    bloch.shift(omega23, laser2_pow)
    bloch.shift(d13real, laser1_det)
    bloch.shift(d23real, laser2_det)
    bloch.shift(d12real, laser12_det)

    results = zeros((9, len(t)))
    vt = v0

    for i in xrange(len(t)):        
        vt, e, hump = spexpv(t[i], bloch, v0)
        vt = matrix(vt).T
        for k in range(0, 9):
            results[k, i] = vt[k]
    return results

if __name__ == "__main__":

    ## Setting up parameters
    # scaling factor to get dimensionless time
    Gt = 2*pi



    # Figure 4.3, page 53
    pl.figure()
    G = ([Gt*1, Gt*0], [Gt*0.99, Gt*0.01], [Gt*0.9, Gt*0.10], [Gt*0.5, Gt*0.5])
    for parami, Gi in enumerate(G):
        # Generate the Bloch-equation matrices of different contributions
        matrices = blochlambda(Gi)
        mnames = ['bloch_atomic', 'laser1_pow', 'laser1_det', 'laser2_pow', 'laser2_det', 'laser12_det']
        for i in range(len(mnames)):
            matrices[i].export_mtx(mnames[i])

        ## Starting parameters
        v0 = matrix(zeros((9, 1)))
        v0[0] = 1

        det13 = 0
        det23 = 0
        omega13 = 1 * Gt
        omega23 = 0 * Gt
        t = linspace(0, 20/Gt, 201)
        args = (mnames, t, v0, det13, det23, omega13, omega23)
        res = timedomain(args)
        t = t * Gt

        pl.subplot(2, 2, parami+1)
        pl.plot(t, res[0, :]*100, 'k-', label="|1>", linewidth=2)
        pl.plot(t, res[1, :]*100, 'k--', label="|2>", linewidth=2)
        pl.plot(t, res[2, :]*100, 'k:', label="|3>", linewidth=2)
        pl.xlabel("Time (1/G)")
        pl.ylabel("State population (%)")
        pl.legend()

    # Figure 4.5, page 55
    pl.figure()
    om23list = [1 * Gt, 2 * Gt, 4 * Gt, 8 * Gt ]
    for parami, omega23 in enumerate(om23list):
        G = [0.9 * Gt, 0.1 * Gt]
        # Generate the Bloch-equation matrices of different contributions
        matrices = blochlambda(G)
        mnames = ['bloch_atomic', 'laser1_pow', 'laser1_det', 'laser2_pow', 'laser2_det', 'laser12_det']
        for i in range(len(mnames)):
            matrices[i].export_mtx(mnames[i])

        ## Starting parameters
        v0 = matrix(zeros((9, 1)))
        v0[0] = 1

        det13 = 0
        det23 = 0
        omega13 = 1 * Gt
        t = linspace(0, 20/Gt, 201)
        args = (mnames, t, v0, det13, det23, omega13, omega23)
        res = timedomain(args)
        t = t * Gt

        pl.subplot(2, 2, parami+1)
        pl.plot(t, res[0, :]*100, 'k-', label="|1>", linewidth=2)
        pl.plot(t, res[1, :]*100, 'k--', label="|2>", linewidth=2)
        pl.plot(t, res[2, :]*100, 'k:', label="|3>", linewidth=2)
        pl.xlabel("Time (1/G)")
        pl.ylabel("State population (%)")
        pl.legend()

    # Figure 4.6, page 56
    pl.figure()
    om23list = [1 * Gt, 2 * Gt, 4 * Gt, 8 * Gt ]
    for parami, omega23 in enumerate(om23list):
        G = [0.9 * Gt, 0.1 * Gt]
        # Generate the Bloch-equation matrices of different contributions
        matrices = blochlambda(G)
        mnames = ['bloch_atomic', 'laser1_pow', 'laser1_det', 'laser2_pow', 'laser2_det', 'laser12_det']
        for i in range(len(mnames)):
            matrices[i].export_mtx(mnames[i])

        ## Starting parameters
        v0 = matrix(zeros((9, 1)))
        v0[0] = 1

        det23 = 0
        omega13 = 1 * Gt
        t = [0, 100 / Gt]
        det13list = linspace(-10*Gt, 10*Gt, 401)
        reslist  = array([])
        for det13 in det13list:
            args = (mnames, t, v0, det13, det23, omega13, omega23)
            res = timedomain(args)
            reslist = append(reslist, 100*res[2, -1])

        pl.subplot(2, 2, parami+1)
        pl.plot(det13list / Gt, reslist, 'k-', label="|3>", linewidth=2)
        pl.ylim([0, max(reslist)*1.05])
        pl.xlim([det13list[0] / Gt, det13list[-1] / Gt])
        pl.xlabel("Time (1/G)")
        pl.ylabel("|3> state population (%)")

    pl.show()
