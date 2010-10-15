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
from numpy import *
from scipy import *
from scipy.linalg import expm
import pylab as pl
from scipy.integrate import ode, quad
from time import time, clock
from pysparse.spmatrix import *
from expv import spexpv

def genmat(G):
    """
    genmat(G)
    Generate the parts of the Bloch-equations to use in the
    actual calculation

    Input:
    G = [G31, G32] the upper level decay rate to the two ground states

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
    

# def timedomain(mnames, split12, det13, tau, reptime, maxpow, npulse, upperlife=1, doplot=False):
def timedomain(args):
    """
    Input:
    mnames, det13, det23, gausstau, gausstheta, gausssigma, trep, npulse, upperlife, doplot

    pulsed laser: laser1
    CW laser: laser2
    """
    mnames, v0, det13, det23, omega23, gausstau, gausstheta, gausssigma, trep, npulse, upperlife, doplot = args
    # mnames, split12, det13, tau, reptime, maxpow, npulse, upperlife, doplot = args
    for name in mnames:
        globals()[name] = ll_mat_from_mtx(name)
    # bloch_atomic, laser1_pow, laser1_det, laser2_pow, laser2_det, laser12_det = matrices
    # detunings
    d13real = det13
    d23real = det23
    d12real = d13real - d23real

    # Between pulses: CW laser is on
    pulseoff = bloch_atomic.copy()
    pulseoff.shift(omega23, laser2_pow)
    pulseoff.shift(d13real, laser1_det)
    pulseoff.shift(d23real, laser2_det)
    pulseoff.shift(d12real, laser12_det)
    pulseofftime = trep - gausstau

    # During pulse - base before adding pulsed power
    # with .shift(omega13(t), laser1_pow)
    pulsebase = pulseoff.copy()
    pulsestep = linspace(0, gausstau, 51)[0:-1]
    pulsedt = pulsestep[1] - pulsestep[0]
    pulsestep = (pulsestep + pulsedt / 2) - gausstau / 2

    # v = matrix(zeros((9, 1)))
    # v[0] = 1/2
    # v[1] = 1/2

    # print pulseoff
    
    results = zeros((9, npulse+1))
    vt = v0
    for k in xrange(0, 9):
        results[k, 0] = vt[k]

    for i in xrange(npulse):

        # pulse on
        v3in = vt[2]
        for pt in pulsestep:
            omega13 = gausspulse(pt, gausstheta, gausssigma)
            pulseon = pulsebase.copy()
            pulseon.shift(omega13, laser1_pow)
            vt, e, hump = spexpv(pulsedt, pulseon, vt)
            vt = matrix(vt).T
        v3diff = vt[2] - v3in
        # print v3diff
        # pulse off
        vt, e, hump = spexpv(pulseofftime, pulseoff, vt)
        vt = matrix(vt).T
        for k in range(0, 9):
            results[k, i+1] = vt[k]

    # if doplot:
    #     tsteps = array(range(0, npulse+1))*reptime

    #     pl.figure()
    #     for i in xrange(3):
    #         pl.plot(tsteps, results[i, :], '-', label="Pop%d"%(i+1))
    #     pl.legend()
    #     pl.title('Population')

    #     pl.figure()
    #     for i in xrange(3):
    #         pl.plot(tsteps, results[i+3, :], '-', label="ChR%d"%(i+1))
    #     pl.legend()
    #     pl.title('Coherence - Real')

    #     pl.figure()
    #     for i in xrange(3):
    #         pl.plot(tsteps, results[i+6, :], '-', label="ChI%d"%(i+1))
    #     pl.legend()
    #     pl.title('Coherence - Imaginary')
    #     pl.show()

    # return results[0:9, -1]
    return results

try:
    import multiprocessing as processing
except:
    import processing

#
# Function run by worker processes
#

def worker(input, output):
    for func, args in iter(input.get, 'STOP'):
        result = calculate(func, args)
        output.put(result)

#
# Function used to calculate result
#

def calculate(func, args):
    res = func(*args)
    return res[2]


# def num_processors():
#     if os.name == 'nt': # Windows
#         return int(os.getenv('NUMBER_OF_PROCESSORS'))
#     else: # glibc (Linux, *BSD, Apple)
#         get_nprocs = ctypes.cdll.libc.get_nprocs
#         get_nprocs.restype = ctypes.c_int
#         get_nprocs.argtypes = []
#         return get_nprocs()

# def pop3(out, trepin, matrices, split12, det13, tau, maxpow, npulse, upperlife):
#     trep = trepin.get()
#     res = timedomain(matrices, split12, det13, tau, trep/split12, maxpow, npulse, upperlife)[2]
#     out.put((trep, res))

# def f(rank, outp, inp, a):
#     x = inp.get()
#     outp.put((x, a*x**2))
#     print rank

if __name__ == "__main__":

    ## Setting up parameters
    # scaling factor to get dimensionless time
    papertime = 1/5.7e6
    Gt = 1
    G = [Gt*0.5, Gt*0.5]

    # Gaussian pulsewidth in dimensionless units
    gausssigma = 100e-15 / papertime
    gausstheta = pi/100
    gausstau = 700e-15 / papertime
    # Repetition rate
    trep = 10e-9 / papertime

    # Generate the Bloch-equation matrices of different contributions
    matrices = genmat(G)
    mnames = ['bloch_atomic', 'laser1_pow', 'laser1_det', 'laser2_pow', 'laser2_det', 'laser12_det']
    for i in range(len(mnames)):
        matrices[i].export_mtx(mnames[i])


    # t = linspace(0, 700, 1401)
    # t0 = 350
    # theta = pi/100
    # sigma = 100

    # print quad(gausspulse, t[0]-t0, t[-1]-t0, args=(theta, sigma))[0]*100
    # pl.plot(t, gausspulse(t-t0, theta, sigma))
    # pl.show()

    ## Starting parameters
    v0 = matrix(zeros((9, 1)))
    v0[0] = 1

    det13 = 0 * Gt
    det23 = 2.5 * Gt
    omega23 = 5 * Gt * 0
    npulse = 100
    args = (mnames, v0, det13, det23, omega23, gausstau, gausstheta, gausssigma, trep, npulse, Gt, False)
    # nump = array(range(0, npulse+1)) / (Gt / (2*pi*trep))
    nump = array(range(0, npulse+1))
    res = timedomain(args)
    pl.figure()
    pl.plot(nump, res[0, :]*100, 's', label="|1>")
    pl.plot(nump, res[1, :]*100, 'o', label="|2>")
    pl.plot(nump, res[2, :]*100, 'x', label="|3>")
    pl.xlabel("Number of pulses")
    pl.ylabel("State population (%)")
    pl.legend()
    pl.figure()
    pl.plot(nump, res[2, :]*100, 'x', label="|3>")
    pl.xlabel("Number of pulses")
    pl.ylabel("Upper state population (%)")
    pl.legend()
    pl.show()

    ########################  Good
    # Gt = 1
    # G31 = Gt * 0.5
    # G32 = Gt - G31
    # Gtot = G31 + G32

    # matrices = genmat([G31, G32])
    # split12 = 1825.9
    # det13 = -split12/4
    # tau = 5e-6
    # maxpow = 300
    # npulse = 1500

    # maxpowlist = [10, 20, 40, 80, 160, 320, 640, 1280]
    # series = 5

 
    
    # treplist = linspace(1.5*2*pi, 3.5*2*pi, 301)
    # NUMBER_OF_PROCESSES = processing.cpu_count()

    # start = time()
    # pool = processing.Pool(processes=NUMBER_OF_PROCESSES)  
    # # result = pool.apply_async(timedomain, (mnames, split12, det13, tau, 6*pi/split12,  maxpow, npulse, Gtot, False))    

    # n = 0
    # for maxpow in maxpowlist:
    #     print "%d / %d : %d" %(n+1, len(maxpowlist), maxpow)
    #     n += 1
    #     TASKS = [(mnames, split12, det13, tau, trep/split12,  maxpow, npulse, Gtot, False) for trep in treplist]
    #     out = pool.map(timedomain, TASKS)
    #     results = zeros((len(treplist), 10))
    #     for i, res in enumerate(out):
    #         results[i, 0] = treplist[i]/2/pi
    #         results[i, 1:10] = res[0:9]
    #     filename = "powercalc%02d_%04d" %(series, maxpow)
    #     savetxt("%s.txt" %filename, results)
    #     for i in xrange(1, 4):
    #         pl.plot(results[:, 0], results[:, i])
    #     pl.savefig("%s.png" %filename)
    #     pl.clf()
    # fin = time() - start
    # print "Totaltime: %f" %fin
    # pool.terminate()
    ########################  End: Good



    # print results
    # for i in xrange(1, 4):
    #     pl.plot(results[:, 0], results[:, i])
    # pl.show()


    # print split12, det13, 10/split12
    # series = True
    # series = False
    # fignum = 2


    # if series:
    #     reptime = 6.5*pi/split12
    #     pops = timedomain(mnames, split12, det13, tau, reptime, maxpow, npulse, upperlife=Gtot, doplot=True)
    # else:
    #     treplist = linspace(3.5*pi, 4.5*pi, 201)
    #     pop1 = array([])
    #     pop2 = array([])
    #     pop3 = array([])
    #     coh1 = array([])
    #     coh2 = array([])
    #     coh3 = array([])
    #     cohi1 = array([])
    #     cohi2 = array([])
    #     cohi3 = array([])
    #     for trep in treplist:
    #         print trep
    #         pops = timedomain(matrices, split12, det13, tau, trep/split12, maxpow, npulse, upperlife=Gtot)
    #         pop1 = append(pop1, pops[0])
    #         pop2 = append(pop2, pops[1])
    #         pop3 = append(pop3, pops[2])
    #         coh1 = append(coh1, pops[3])
    #         coh2 = append(coh2, pops[4])
    #         coh3 = append(coh3, pops[5])
    #         cohi1 = append(cohi1, pops[6])
    #         cohi2 = append(cohi2, pops[7])
    #         cohi3 = append(cohi3, pops[8])
    #     treplist = treplist/2/pi
    #     pl.figure()
    #     # pl.title("Population")
    #     pl.plot(treplist, pop1, label="Pop |1>")
    #     pl.plot(treplist, pop2, label="Pop |2>")
    #     pl.plot(treplist, pop3, label="Pop |3>")
    #     pl.xlabel('Repetition rate divider')
    #     pl.ylabel('Population')
    #     pl.legend()
    #     pl.savefig("reprate%02d_1.pdf" %fignum)
    #     pl.figure()
    #     pl.title("Real coherence")
    #     pl.plot(treplist, coh1, label="1")
    #     pl.plot(treplist, coh2, label="2")
    #     pl.plot(treplist, coh3, label="3")
    #     pl.legend()
    #     pl.savefig("reprate%02d_2.pdf" %fignum)
    #     pl.figure()
    #     pl.title("Imaginary coherence")
    #     pl.plot(treplist, cohi1, label="1")
    #     pl.plot(treplist, cohi2, label="2")
    #     pl.plot(treplist, cohi3, label="3")
    #     pl.legend()
    #     pl.savefig("reprate%02d_3.pdf" %fignum)
    #     pl.show()
