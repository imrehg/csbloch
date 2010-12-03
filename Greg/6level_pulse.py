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
    mnames, v, t, omega, laser_det, delta = args
    for name in mnames:
        globals()[name] = ll_mat_from_mtx(name)
    # bloch_atomic, laser_pow = matrices

    # v = matrix(zeros((36, 1)))
    # v[1] = 1

    laseron = bloch_atomic.copy()
    laseron.shift(omega, laser_pow)
    laseron.shift(1, laser_det(delta))
    # print laseron
    # print laseron

    t = np.linspace(0, t, 1001)
    results = zeros((len(v), len(t)))
    vt = v
    for k in xrange(0, 36):
        results[k, 0] = vt[k]
    for i in xrange(1, len(t)):
        vt, e, hump = spexpv(t[i]-t[i-1], laseron, vt)
        vt = matrix(vt).T

        for k in range(0, 36):
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

    print "Reptime: %f, detuning: %f" %(trep, delta)

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

    TASKS = [(kargs['mnames'], kargs['v'], bloch.laser_detd2,
              kargs['detu'], kargs['gausspulse'], kargs['gausssigma'],
              kargs['gausstheta'], kargs['gausstau'], 30, trep, kargs['npulse']) for trep in treplist]
    try:
        out = pool.map_async(pulsescalc, TASKS).get(999999)
    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(0)
    nelems = len(kargs['v'])
    results = zeros((len(treplist), nelems+1))
    for i, res in enumerate(out):
            results[i, 0] = treplist[i]
            for k in xrange(nelems):
                results[i, k+1] = res[k][-1]

    d01 = abs(bloch.laser_detd2(0)[21,6])
    timescale = d01/2/np.pi
    pl.figure()
    for i in xrange(2,6):
        pl.plot(results[:, 0]*timescale, results[:, i+1], '.-', label="|%d>"%i)
    pl.plot(results[:, 0]*timescale, results[:, 3]+results[:, 4]+results[:, 5]+results[:, 6], 'k--', label='total')
    pl.legend(loc='best')
    pl.title('Upper state populations')
    pl.xlim([results[0,0]*timescale, results[-1,0]*timescale])
    pl.xlabel('Repetition rate divider')
    pl.ylabel('Population')
    if savefigs:
        figname = "%s_pop.png" %savefigs
        pl.savefig(figname)
        # catchimgup("Population", "%s" %figname)

    pl.figure()
    for i in [6]:
        pl.plot(results[:,0]*timescale, results[:, i+1], '--', label='Real')
    for i in [21]:
        pl.plot(results[:,0]*timescale, results[:, i+1], ':', label='Imag')
    pl.legend(loc='best')
    pl.title('Ground state coherence')
    pl.xlim([results[0,0]*timescale, results[-1,0]*timescale])
    pl.xlabel('Repetition rate divider')
    pl.ylabel('Coherence')
    if savefigs:
        figname = "%s_coh.png" %savefigs
        pl.savefig(figname)
        # catchimgup("Coherence", "%s" %figname)

def detuscan(detulist, pool, savefigs=None, **kargs):


    # # mnames, v, laser_det, delta, pulsefunc, pulseshape, pulsearea, pulsetotaltime, pulsepoints, trep, npulse = args

    TASKS = [(kargs['mnames'], kargs['v'], bloch.laser_detd2,
              detu, kargs['gausspulse'], kargs['gausssigma'],
              kargs['gausstheta'], kargs['gausstau'], 30, kargs['trep'], kargs['npulse']) for detu in detulist]
    try:
        out = pool.map_async(pulsescalc, TASKS).get(999999)
    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(0)
    nelems = len(kargs['v'])
    results = zeros((len(detulist), nelems+1))
    for i, res in enumerate(out):
            results[i, 0] = detulist[i]
            for k in xrange(nelems):
                results[i, k+1] = res[k][-1]
    if savefigs:
        filename = "%s_results.txt" %savefigs
        savetxt(filename, results)

    pl.figure()
    for i in xrange(2,6):
        pl.plot(results[:, 0], results[:, i+1], '.-', label="|%d>"%i)
    pl.plot(results[:, 0], results[:, 3]+results[:, 4]+results[:, 5]+results[:, 6], 'k--', label='total')
    pl.legend(loc='best')
    pl.title('Upper state populations')
    pl.xlim([results[0,0], results[-1,0]])
    if savefigs:
        figname = "%s_pop.png" %savefigs
        pl.savefig(figname)
        # catchimgup("Population", "%s" %figname)

    pl.figure()
    for i in [6]:
        pl.plot(results[:,0], results[:, i+1], '--', label='Real')
    for i in [21]:
        pl.plot(results[:,0], results[:, i+1], ':', label='Imag')
    pl.legend(loc='best')
    pl.title('Ground state coherence')
    pl.xlim([results[0,0], results[-1,0]          ])
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

    # bloch.laser_detd2(0)
    # sys.exit(0)

    matrices = bloch.blochd2()
    mnames = ['bloch_atomic', 'laser_pow']
    for i in range(len(mnames)):
        matrices[i].export_mtx(mnames[i])

    v = matrix(zeros((36, 1)))
    # v[0] = 0.5
    # v[1] = 0.5
    v[2] = 1

    
    # TASKS = [(mnames, t, bloch.laser_det) for t in range(10,1000)]
    # out = pool.map(timedomain, TASKS)
    # print out[0]
    # print bloch.laser_det(255)

    # omega = np.pi
    # delta = -1756
    # maxt = 10
    # params = [mnames, v, maxt, omega, bloch.laser_detd2, delta]
    # res = timedomain(params)
    # t = linspace(0, maxt, 1001)
    # # print res[2, :]
    # # markers = ['k-', 'b--', 'r-.', 'g:']
    # # for i in range(3, -1, -1):
    # #     pl.plot(t, res[i, :], markers[i], label='|%d>' % (i), linewidth=2)
    # [pl.plot(t, res[i, :]) for i in range(6)]
    # # pl.legend(loc='best')
    # pl.show()
    # sys.exit(0)



    # #### D2 line frequency comb - time dependent
    # G = 32.889e6
    # timescale = 1/G
    # # pulse length = 1ps
    # tauP = 1e-12 / timescale
    # gausssigma = tauP * sqrt(pi/(2*log(2)))
    # gausstheta = 1/50
    # gausstau = 4 * gausssigma

    # v = matrix(zeros((36, 1)))
    # v[0] = 0.5
    # v[1] = 0.5

    # # 
    # d01 = abs(bloch.laser_detd2(0)[21,6])
    # trep = 100.5/(d01/2/np.pi)
    # # print("trep %f" %trep)
    # # print("gtau %f" %gausstau)
    # # sys.exit(0)
    # npulse = 100
    # # detu = -2009.3/2
    # detu = 0

    # # mnames, v, laser_det, delta, pulsefunc, pulseshape, pulsearea, pulsetotaltime, pulsepoints, trep, npulse = args
    # res =  pulsescalc((mnames, v, bloch.laser_detd2, detu, gausspulse, gausssigma, gausstheta, gausstau, 30, trep, npulse))
    # pulses = range(npulse+1)
    # pl.figure(1)
    # for i in xrange(0,2):
    #     pl.plot(pulses, res[i, :], '.-', label='%d'%i)
    # pl.plot(pulses, res[0, :]+res[1, :], 'kx-', label='Sum')
    # pl.title('Lower state populations')
    # pl.legend(loc='best')

    # pl.figure(0)
    # for i in xrange(2,6):
    #     pl.plot(pulses, res[i, :], '.-', label='%d'%i)
    # pl.plot(pulses, res[2, :]+res[3, :]+res[4, :]+res[5, :], 'kx-', label='Sum')
    # pl.title('Upper state populations')
    # pl.legend(loc='best')

    # # pl.figure()
    # # for i in xrange(4, 6):
    # #     pl.plot(pulses, res[i, :], label='%d'%i)
    # # pl.title('Real Coherence')
    # # pl.legend(loc='best')
    # # pl.figure()
    # # for i in xrange(10, ):
    # #     pl.plot(pulses, res[i, :], label='%d'%i)
    # # pl.title('Imag Coherence')
    # # pl.legend(loc='best')
    # pl.show()
    # sys.exit(0)



    # sys.exit(0)
    # # #### Fixed repetition rate and scanned detuning
    G = 32.889e6
    timescale = 1/G
    # pulse length = 1ps
    tauP = 1e-12 / timescale
    gausssigma = tauP * sqrt(pi/(2*log(2)))
    gausstheta = 1/25
    gausstau = 4 * gausssigma

    v = matrix(zeros((36, 1)))
    v[0] = 0.5
    v[1] = 0.5

    d01 = abs(bloch.laser_detd2(0)[21,6])
    # trep = 100.5/(d01/2/np.pi)
    # print("trep %f" %trep)
    # print("gtau %f" %gausstau)
    # sys.exit(0)
    npulse = 500
    # # detu = -2009.3/2
    # detu = 0

    trepbase = 100
    dt = 0.04
    ntrep = 81

    treplist = linspace(trepbase-dt, trepbase+dt, ntrep)/(d01/2/np.pi)
    figs = ['6level_f32',
            '6level_f33',
            '6level_f34',
            '6level_f35',
            '6level_f42',
            '6level_f43',
            '6level_f44',
            '6level_f45']
    detus = [0, 28.89, 28.89+38.46, 28.89+38.46+47.97,
             -1756.33, -1756.33+28.89, -1756.33+28.89+38.46, -1756.33+28.89+38.46+47.97]

    # loop through all the settings
    for pars in zip(figs, detus):
        repratescan(treplist, pool=pool, savefigs=pars[0], mnames=mnames, v=v,
                detu=pars[1], gausspulse=gausspulse, gausssigma=gausssigma,
                gausstheta=gausstheta, gausstau=gausstau,
                npulse=npulse)


    # # on resonance with F=4->F'5
    # detu = -1638.01
    # repratescan(treplist, pool=pool, savefigs="6level_f45", mnames=mnames, v=v,
    #             detu=detu, gausspulse=gausspulse, gausssigma=gausssigma,
    #             gausstheta=gausstheta, gausstau=gausstau,
    #             npulse=npulse)
    pl.show()

    # trep = 1/80e6/timescale
    # detulist = 2*pi*linspace(-1, 1, ndetu)*80e6/G
    # detuscan(detulist, pool=pool, savefigs="4level_detu_f3", mnames=mnames, v=v,
    #             gausspulse=gausspulse, gausssigma=gausssigma,
    #             gausstheta=gausstheta, gausstau=gausstau, trep=trep,
    #             npulse=npulse)
    # pl.show()

    # # #### Fixed repetition rate and scanned detuning
    # # print matrices[1]
    # # sys.exit(0)
    # d01 = abs(bloch.laser_det(0)[10,4])

    # G = 28.743e6
    # timescale = 1/G
    # tauP = 1e-12 / timescale
    # gausssigma = tauP * sqrt(pi/(2*log(2)))
    # gausstheta = 1/100
    # gausstau = 4 * gausssigma


    # detu = 0
    # # trep = 1/80e6/timescale
    # # trep = 1/G/timescale/100
    # trep = 80.5/(d01/2/np.pi)
    # # print("trep %f" %trep)
    # # print("gtau %f" %gausstau)
    # npulse = 1200

    # dt = 0.1
    # ntrep = 151
    # ndetu = 301

    # trep = 1/80e6/timescale
    # detulist = 2*pi*linspace(-1, 1, ndetu)*80e6/G
    # detuscan(detulist, pool=pool, savefigs="4level_detu_f3", mnames=mnames, v=v,
    #             gausspulse=gausspulse, gausssigma=gausssigma,
    #             gausstheta=gausstheta, gausstau=gausstau, trep=trep,
    #             npulse=npulse)
    # pl.show()


    # ##### Timedomain
    # pulseargs = (mnames, v, bloch.laser_det, detu, gausspulse, gausssigma, gausstheta, gausstau, 30, trep, npulse)
    # res = pulsescalc(pulseargs)
    # for i in [0, 1, 2, 3]:
    #     pl.plot(res[i, :])
    # pl.plot(res[0,:]+res[1,:]+res[2,:]+res[3,:])
    # pl.plot(res[2,:]+res[3,:], 'k--')
    # pl.show()
    # sys.exit(0)


    # ##### Repratespectro
    # treplist = linspace(80-dt, 80+dt, ntrep)/(d01/2/np.pi)
    # # on resonance with F=3->F'4
    # detu = 0
    # repratescan(treplist, pool=pool, savefigs="4level_f3", mnames=mnames, v=v,
    #             detu=detu, gausspulse=gausspulse, gausssigma=gausssigma,
    #             gausstheta=gausstheta, gausstau=gausstau,
    #             npulse=npulse)

    # # halfway between F=3,4
    # detu = -d01/2
    # repratescan(treplist, pool=pool, savefig="4level_f34", mnames=mnames, v=v,
    #             detu=detu, gausspulse=gausspulse, gausssigma=gausssigma,
    #             gausstheta=gausstheta, gausstau=gausstau,
    #             npulse=npulse)

    # # on resonance with F=4->F'=3
    # detu = -d01
    # repratescan(treplist, pool=pool, savefig="4level_f4", mnames=mnames, v=v,
    #             detu=detu, gausspulse=gausspulse, gausssigma=gausssigma,
    #             gausstheta=gausstheta, gausstau=gausstau,
    #             npulse=npulse)

    # # Detuning with half-gap with respect to the repetition rate
    # detu = -d01/80/2
    # repratescan(treplist, pool=pool, savefig="4level_half", mnames=mnames, v=v,
    #             detu=detu, gausspulse=gausspulse, gausssigma=gausssigma,
    #             gausstheta=gausstheta, gausstau=gausstau,
    #             npulse=npulse)
    # pl.show()

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
