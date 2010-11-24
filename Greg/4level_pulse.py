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
from catchupload import *

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
    v[1] = 1

    laseron = bloch_atomic.copy()
    laseron.shift(omega, laser_pow)
    laseron.shift(1, laser_det(delta))
    # print laseron
    # print laseron

    t = np.linspace(0, t, 1001)
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

    print trep

    # Between pulses: CW laser is off
    pulseoff = bloch_atomic.copy()
    pulseoff.shift(1, laser_det(delta))
    pulseofftime = trep - pulsetotaltime


    # During pulse - base before adding pulsed power
    pulsebase = pulseoff.copy()

    pulsestep = linspace(0, pulsetotaltime, pulsepoints)[0:-1]
    pulsedt = pulsestep[1] - pulsestep[0]
    pulsestep = (pulsestep + pulsedt / 2) - pulsetotaltime / 2
    pulseblochs = []
    omegas = []
    for pt in pulsestep:
        omega = pulsefunc(pt, pulsearea, pulseshape)
        pulseon = pulsebase.copy()
        pulseon.shift(omega, laser_pow)
        pulseblochs += [pulseon]
        omegas += [omega]

    nelems = len(v0)
    results = zeros((nelems, npulse+1))
    vt = v0
    for k in xrange(0, nelems):
        results[k, 0] = vt[k]

    for i in xrange(npulse):
        # pulse on
        # v3in = vt[2]
        # for p in xrange(len(pulsestep)):
        for nextstep in pulseblochs:
            vt, e, hump = spexpv(pulsedt, nextstep, vt)
            vt = matrix(vt).T
        # pulse off
        vt, e, hump = spexpv(pulseofftime, pulseoff, vt)
        vt = matrix(vt).T
        for k in xrange(0, nelems):
            results[k, i+1] = vt[k]
    return results



def repratescan(treplist, pool, savefigs=None, **kargs):


    # # mnames, v, laser_det, delta, pulsefunc, pulseshape, pulsearea, pulsetotaltime, pulsepoints, trep, npulse = args

    TASKS = [(kargs['mnames'], kargs['v'], bloch.laser_det,
              kargs['detu'], kargs['gausspulse'], kargs['gausssigma'],
              kargs['gausstheta'], kargs['gausstau'], 30, trep, kargs['npulse']) for trep in treplist]
    try:
        out = pool.map_async(pulsescalc, TASKS).get(1000)
    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(0)
    results = zeros((len(treplist), 17))
    for i, res in enumerate(out):
            results[i, 0] = treplist[i]
            for k in xrange(16):
                results[i, k+1] = res[k][-1]

    timescale = 2009.3/2/np.pi
    pl.figure()
    for i in xrange(2,4):
        pl.plot(results[:, 0]*timescale, results[:, i+1], '.-', label="|%d>"%i)
    pl.plot(results[:, 0]*timescale, results[:, 3]+results[:, 4], 'k--', label='total')
    pl.legend(loc='best')
    pl.title('Upper state populations')
    if savefigs:
        figname = "%s_pop.png" %savefigs
        pl.savefig(figname)
        # catchimgup("Population", "%s" %figname)

    pl.figure()
    for i in [4]:
        pl.plot(results[:,0]*timescale, results[:, i+1], '--', label='Real')
    for i in [10]:
        pl.plot(results[:,0]*timescale, results[:, i+1], ':', label='Imag')
    pl.legend(loc='best')
    pl.title('Ground state coherence')
    if savefigs:
        figname = "%s_coh.png" %savefigs
        pl.savefig(figname)
        # catchimgup("Coherence", "%s" %figname)


def main():

    try:
        import multiprocessing as processing
    except:
        import processing
    NUMBER_OF_PROCESSES = processing.cpu_count()
    pool = processing.Pool(processes=NUMBER_OF_PROCESSES)  


    matrices = bloch.blochd1()
    # print matrices[1]
    # print bloch.laser_det(0)
    # print bloch.laser_det(10)
    # sys.exit(0)
    mnames = ['bloch_atomic', 'laser_pow']
    for i in range(len(mnames)):
        matrices[i].export_mtx(mnames[i])

    # print matrices[0]
    # print matrices[1]
    # sys.exit(0)
    v = matrix(zeros((16, 1)))
    v[0] = 0.5
    v[1] = 0.5

    
    # TASKS = [(mnames, t, bloch.laser_det) for t in range(10,1000)]
    # out = pool.map(timedomain, TASKS)
    # print out[0]
    # print bloch.laser_det(255)

    # omega = 2*pi
    # delta = -1755
    # print bloch.laser_det(delta)
    # maxt = 10
    # params = [mnames, maxt, omega, bloch.laser_det, delta]
    # res = timedomain(params)
    # t = linspace(0, 10, 1001)
    # # print res[2, :]
    # markers = ['k-', 'b--', 'r-.', 'g:']
    # for i in range(3, -1, -1):
    #     pl.plot(t, res[i, :], markers[i], label='|%d>' % (i), linewidth=2)
    # pl.legend(loc='best')
    # pl.show()

    # #### new
    # G = 28.743e6
    # timescale = 1/G
    # tauP = 1e-12 / timescale
    # gausssigma = tauP * sqrt(pi/(2*log(2)))
    # gausstheta = 1/50
    # gausstau = 4 * gausssigma

    # trep = 1/abs(bloch.laser_det(0)[10,4])
    # print("trep %f" %trep)
    # print("gtau %f" %gausstau)
    # npulse = 500
    # # detu = -2009.3/2
    # detu = 0

    # # mnames, v, laser_det, delta, pulsefunc, pulseshape, pulsearea, pulsetotaltime, pulsepoints, trep, npulse = args
    # res =  pulsescalc((mnames, v, bloch.laser_det, detu, gausspulse, gausssigma, gausstheta, gausstau, 30, trep, npulse))
    # pulses = range(npulse+1)
    # pl.figure(1)
    # for i in xrange(0,2):
    #     pl.plot(pulses, res[i, :], '.-', label='%d'%i)
    # pl.plot(pulses, res[0, :]+res[1, :], 'kx-', label='Sum')
    # pl.title('Lower state populations')
    # pl.legend(loc='best')

    # pl.figure(0)
    # for i in xrange(2,4):
    #     pl.plot(pulses, res[i, :], '.-', label='%d'%i)
    # pl.plot(pulses, res[2, :]+res[3, :], 'kx-', label='Sum')
    # pl.title('Upper state populations')
    # pl.legend(loc='best')

    # pl.figure()
    # for i in xrange(4, 5):
    #     pl.plot(pulses, res[i, :], label='%d'%i)
    # pl.title('Real Coherence')
    # pl.legend(loc='best')
    # pl.figure()
    # for i in xrange(10, 11):
    #     pl.plot(pulses, res[i, :], label='%d'%i)
    # pl.title('Imag Coherence')
    # pl.legend(loc='best')
    # pl.show()


    # #### new22
    # print matrices[1]
    # sys.exit(0)
    d01 = abs(bloch.laser_det(0)[10,4])

    G = 28.743e6
    timescale = 1/G
    tauP = 1e-12 / timescale
    gausssigma = tauP * sqrt(pi/(2*log(2)))
    gausstheta = np.pi/10
    gausstau = 4 * gausssigma


    detu = 0
    # trep = 1/80e6/timescale
    # trep = 1/G/timescale/100
    trep = 80.5/(d01/2/np.pi)
    # print("trep %f" %trep)
    # print("gtau %f" %gausstau)
    npulse = 100

    dt = 0.5
    ntrep = 301

    # pulseargs = (mnames, v, bloch.laser_det, detu, gausspulse, gausssigma, gausstheta, gausstau, 30, trep, npulse)
    # res = pulsescalc(pulseargs)
    # for i in [0, 1, 2, 3]:
    #     pl.plot(res[i, :])
    # pl.plot(res[0,:]+res[1,:]+res[2,:]+res[3,:])
    # pl.plot(res[2,:]+res[3,:], 'k--')
    # pl.show()
    # sys.exit(0)

    treplist = linspace(80-dt, 80+dt, ntrep)/(d01/2/np.pi)
    detu = 0
    repratescan(treplist, pool=pool, savefigs="4level_f3", mnames=mnames, v=v,
                detu=detu, gausspulse=gausspulse, gausssigma=gausssigma,
                gausstheta=gausstheta, gausstau=gausstau,
                npulse=npulse)
    detu = -d01/2
    repratescan(treplist, pool=pool, savefig="4level_f34", mnames=mnames, v=v,
                detu=detu, gausspulse=gausspulse, gausssigma=gausssigma,
                gausstheta=gausstheta, gausstau=gausstau,
                npulse=npulse)
    detu = -d01
    repratescan(treplist, pool=pool, savefig="4level_f4", mnames=mnames, v=v,
                detu=detu, gausspulse=gausspulse, gausssigma=gausssigma,
                gausstheta=gausstheta, gausstau=gausstau,
                npulse=npulse)
    pl.show()

    # detulist = np.linspace(-300, 0, 41)

    # print (d01+detu)/(2*pi/trep)

    # # mnames, v, laser_det, delta, pulsefunc, pulseshape, pulsearea, pulsetotaltime, pulsepoints, trep, npulse = args
    # res =  pulsescalc((mnames, v, bloch.laser_det, detu, gausspulse, gausssigma, gausstheta, gausstau, 30, trep, npulse))
    # pulses = range(npulse+1)
    # pl.figure(1)
    # for i in xrange(0,2):
    #     pl.plot(pulses, res[i, :], '.-', label='%d'%i)
    # pl.plot(pulses, res[0, :]+res[1, :], 'kx-', label='Sum')
    # pl.title('Lower state populations')
    # pl.legend(loc='best')

    # pl.figure(0)
    # for i in xrange(2,4):
    #     pl.plot(pulses, res[i, :], '.-', label='%d'%i)
    # pl.plot(pulses, res[2, :]+res[3, :], 'kx-', label='Sum')
    # pl.title('Upper state populations')
    # pl.legend(loc='best')

    # pl.figure()
    # for i in xrange(4, 5):
    #     pl.plot(pulses, res[i, :], label='%d'%i)
    # pl.title('Real Coherence')
    # pl.legend(loc='best')
    # pl.figure()
    # for i in xrange(10, 11):
    #     pl.plot(pulses, res[i, :], label='%d'%i)
    # pl.title('Imag Coherence')
    # pl.legend(loc='best')
    # pl.show()


    # TASKS = [(mnames, v, bloch.laser_det, detu, gausspulse, gausssigma, gausstheta, gausstau, 30, trep, npulse) for detu in detulist]
    # out = pool.map_async(pulsescalc, TASKS).get(1000)
    # results = zeros((len(detulist), 16))
    # for i, res in enumerate(out):
    #         results[i, 0] = detulist[i]
    #         for k in xrange(9):
    #             results[i, k+1] = res[k][-1]

    # pl.figure(1)
    # for i in xrange(3,5):
    #     pl.plot(results[:, 0], results[:, i], '.-')
    # pl.title('Upper state populations')

    # pl.figure(2)
    # for i in xrange(4, 6):
    #     pl.plot(pulses, res[i, :], '.-')
    # pl.title('Coherence')

    
    # # mnames, v, laser_det, delta, pulsefunc, pulseshape, pulsearea, pulsetotaltime, pulsepoints, trep, npulse = args
    # res =  pulsescalc((mnames, v, bloch.laser_det, detu, gausspulse, gausssigma, gausstheta, gausstau, 30, trep, npulse))
    # pulses = range(npulse+1)
    # pl.figure(1)
    # for i in xrange(0,2):
    #     pl.plot(pulses, res[i, :], '.-', label='%d'%i)
    # pl.plot(pulses, res[0, :]+res[1, :], 'kx-', label='Sum')
    # pl.title('Lower state populations')
    # pl.legend(loc='best')

    # pl.figure(0)
    # for i in xrange(2,4):
    #     pl.plot(pulses, res[i, :], '.-', label='%d'%i)
    # pl.plot(pulses, res[2, :]+res[3, :], 'kx-', label='Sum')
    # pl.title('Upper state populations')
    # pl.legend(loc='best')

    # pl.figure()
    # for i in xrange(4, 5):
    #     pl.plot(pulses, res[i, :], label='%d'%i)
    # pl.title('Real Coherence')
    # pl.legend(loc='best')
    # pl.figure()
    # for i in xrange(10, 11):
    #     pl.plot(pulses, res[i, :], label='%d'%i)
    # pl.title('Imag Coherence')
    # pl.legend(loc='best')
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


    # trep = 1/80e6/timescale

    # npulse = 200
    # # detu = -2009.3/2

    # # mnames, v, laser_det, delta, pulsefunc, pulseshape, pulsearea, pulsetotaltime, pulsepoints, trep, npulse = args
    # detulist = np.linspace(-300e6, 300e6, 31)*2*np.pi/G

    # TASKS = [(mnames, v, bloch.laser_det, detu, gausspulse, gausssigma, gausstheta, gausstau, 21, trep, npulse) for detu in detulist]
    # out = pool.map_async(pulsescalc, TASKS).get(999999)
    # results = zeros((len(detulist), 16))
    # for i, res in enumerate(out):
    #         results[i, 0] = detulist[i]*G/2/np.pi/1e6
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

    # pl.show()

if __name__ == "__main__":
    main()
