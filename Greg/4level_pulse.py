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


def gausspulse(t, theta, sigma):
    """
    Gaussian pulse shape:
    gausspulse(t, theta, sigma)

    Input:
    t = requested times
    theta = pulse area
    sigma = pulse width

    Output:
    as in the paper's Eq. (7)
    """
    return (theta/sigma)*exp(-pi*(t/sigma)**2)


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


def pulsescalc(args):
    """
    Calculate time evolution with pulses
    """
    # mnames, split12, det13, tau, reptime, maxpow, npulse, upperlife, doplot = args
    mnames, v0, laser_det, delta, pulsefunc, pulseshape, pulsearea, pulsetotaltime, pulsepoints, trep, npulse = args
    for name in mnames:
        globals()[name] = ll_mat_from_mtx(name)
    # bloch_atomic, laser_pow = matrices

    # Between pulses: CW laser is off
    pulseoff = bloch_atomic.copy()
    pulseofftime = trep - pulsetotaltime

    # During pulse - base before adding pulsed power
    pulsebase = bloch_atomic.copy()
    pulsebase.shift(1, laser_det(delta))
    pulsestep = linspace(0, pulsetotaltime, pulsepoints)[0:-1]
    pulsedt = pulsestep[1] - pulsestep[0]
    pulsestep = (pulsestep + pulsedt / 2) - pulsetotaltime / 2
    pulseblochs = []
    for pt in pulsestep:
        omega = pulsefunc(pt, pulsearea, pulseshape)
        pulseon = pulsebase.copy()
        pulseon.shift(omega, laser_pow)
        pulseblochs += [pulseon]

    results = zeros((9, npulse+1))
    vt = v0
    for k in xrange(0, 9):
        results[k, 0] = vt[k]

    for i in xrange(npulse):
        # pulse on
        v3in = vt[2]
        for p in xrange(len(pulsestep)):
            vt, e, hump = spexpv(pulsedt, pulseblochs[p], vt)
            vt = matrix(vt).T
        # pulse off
        vt, e, hump = spexpv(pulseofftime, pulseoff, vt)
        vt = matrix(vt).T
        for k in range(0, 9):
            results[k, i+1] = vt[k]
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

    
    v = matrix(zeros((16, 1)))
    v[0] = 0.5
    v[1] = 0.5

    # TASKS = [(mnames, t, bloch.laser_det) for t in range(10,1000)]
    # out = pool.map(timedomain, TASKS)
    # print out[0]
    # print bloch.laser_det(255)

    # omega = 2*pi
    # delta = 0
    # maxt = 10
    # params = [mnames, maxt, omega, bloch.laser_det, delta]
    # res = timedomain(params)
    # t = linspace(0, 10)
    # # print res[2, :]
    # markers = ['k-', 'b--', 'r-.', 'g:']
    # for i in range(3, -1, -1):
    #     pl.plot(t, res[i, :], markers[i], label='|%d>' % (i), linewidth=2)
    # pl.legend(loc='best')
    # pl.show()

    G = 28.743e6
    timescale = 1/G
    tauP = 1e-12 / timescale
    gausssigma = tauP * sqrt(pi/(2*log(2)))
    gausstheta = 1/100
    gausstau = 8 * gausssigma

    # trep = 100/2009.3

    # npulse = 200
    # detu = -2009.3/2
    # # mnames, v, laser_det, delta, pulsefunc, pulseshape, pulsearea, pulsetotaltime, pulsepoints, trep, npulse = args
    # res =  pulsescalc((mnames, v, bloch.laser_det, detu, gausspulse, gausssigma, gausstheta, gausstau, 21, trep, npulse))

    # pulses = range(npulse+1)
    # pl.figure(1)
    # for i in xrange(2,4):
    #     pl.plot(pulses, res[i, :], '.-')
    # pl.plot(pulses, res[2, :]+res[3, :], 'kx-')
    # pl.title('Upper state populations')

    # pl.figure(2)
    # for i in xrange(4, 6):
    #     pl.plot(pulses, res[i, :], '.-')
    # pl.title('Coherence')
    # pl.show()

    # ##### Repetition rate
    # treplist = linspace(95.5, 100.5, 31)/2009.3
    # npulse = 200
    # detu = -2009.3/2

    # # mnames, v, laser_det, delta, pulsefunc, pulseshape, pulsearea, pulsetotaltime, pulsepoints, trep, npulse = args

    # TASKS = [(mnames, v, bloch.laser_det, detu, gausspulse, gausssigma, gausstheta, gausstau, 21, trep, npulse) for trep in treplist]
    # out = pool.map(pulsescalc, TASKS)
    # results = zeros((len(treplist), 16))
    # for i, res in enumerate(out):
    #         results[i, 0] = treplist[i]
    #         for k in xrange(9):
    #             results[i, k+1] = res[k][-1]

    # pl.figure(1)
    # for i in xrange(3,5):
    #     pl.plot(results[:, 0], results[:, i], '.-')
    # pl.title('Upper state populations')

    # # pl.figure(2)
    # # for i in xrange(4, 6):
    # #     pl.plot(pulses, res[i, :], '.-')
    # # pl.title('Coherence')


    trep = 1/80e6/timescale

    npulse = 200
    # detu = -2009.3/2

    # mnames, v, laser_det, delta, pulsefunc, pulseshape, pulsearea, pulsetotaltime, pulsepoints, trep, npulse = args
    detulist = np.linspace(-300e6, 300e6, 31)*2*np.pi/G

    TASKS = [(mnames, v, bloch.laser_det, detu, gausspulse, gausssigma, gausstheta, gausstau, 21, trep, npulse) for detu in detulist]
    out = pool.map_async(pulsescalc, TASKS).get(999999)
    results = zeros((len(detulist), 16))
    for i, res in enumerate(out):
            results[i, 0] = detulist[i]*G/2/np.pi/1e6
            for k in xrange(9):
                results[i, k+1] = res[k][-1]

    pl.figure(1)
    for i in xrange(3,5):
        pl.plot(results[:, 0], results[:, i], '.-')
    pl.title('Upper state populations')

    # pl.figure(2)
    # for i in xrange(4, 6):
    #     pl.plot(pulses, res[i, :], '.-')
    # pl.title('Coherence')

    pl.show()

if __name__ == "__main__":
    main()
