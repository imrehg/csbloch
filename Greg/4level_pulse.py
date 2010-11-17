#!/usr/bin/env python

from __future__ import division
import numpy as np
from scipy import *
from scipy.linalg import expm
import pylab as pl
from scipy.integrate import ode
from pysparse.spmatrix import *
import physicspy.quantum as quantum
import blochmatrix as bloch
from expv import spexpv
import sys

def timedomain(args):
    # mnames, split12, det13, tau, reptime, maxpow, npulse, upperlife, doplot = args
    mnames, t, omega, laser_det, delta = args
    for name in mnames:
        globals()[name] = ll_mat_from_mtx(name)
    # bloch_atomic, laser_pow = matrices

    v = matrix(zeros((16, 1)))
    v[0] = 0.5
    v[1] = 0.5

    laseron = bloch_atomic.copy()
    laseron.shift(omega, laser_pow)
    laseron.shift(1, laser_det(delta))
    # print laseron

    t = np.linspace(0, t)
    results = zeros((16, len(t)))
    vt = v
    for k in xrange(0, 16):
        results[k, 0] = vt[k]
    for i in xrange(1, len(t)):
        vt, e, hump = spexpv(t[i]-t[i-1], laseron, vt)
        vt = matrix(vt).T

        for k in range(0, 16):
            results[k, i] = vt[k]

    return results


def main():

    try:
        import multiprocessing as processing
    except:
        import processing
    NUMBER_OF_PROCESSES = processing.cpu_count()
    pool = processing.Pool(processes=NUMBER_OF_PROCESSES)  


    matrices = bloch.blochd1()
    mnames = ['bloch_atomic', 'laser_pow']
    for i in range(len(mnames)):
        matrices[i].export_mtx(mnames[i])

    

    # TASKS = [(mnames, t, bloch.laser_det) for t in range(10,1000)]
    # out = pool.map(timedomain, TASKS)
    # print out[0]
    # print bloch.laser_det(255)
    omega = 1
    delta = 0
    maxt = 10
    params = [mnames, maxt, omega, bloch.laser_det, delta]
    res = timedomain(params)
    t = linspace(0, 10)
    # print res[2, :]
    markers = ['k-', 'b--', 'r-.', 'g:']
    for i in range(3, -1, -1):
        pl.plot(t, res[i, :], markers[i], label='|%d>' % (i), linewidth=2)
    pl.legend(loc='best')
    pl.show()

if __name__ == "__main__":
    main()
