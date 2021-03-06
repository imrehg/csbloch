#!/usr/bin/env python

from __future__ import division
from numpy import *
from scipy import *
from scipy.linalg import expm
import pylab as pl
from scipy.integrate import ode
from time import time
from pysparse.spmatrix import ll_mat
from expv import spexpv

# Level structure:
#     ---- (3)
#     /  \
# ----(1) ----(2)

# # Upper level lifetime, detunings and linewidths
# Gamma is in anglular units - e.g. 1 MHz => 2pi x 1 MHz
Gt = 1
d13 = 0
G31 = Gt * 0.9
d23 = 0
G32 = Gt - G31
Gtot = G31 + G32

om13real = 0.5
om23real = 2
d13real = -1.3
d23real = -1.3
start = 0

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


# print out
# nt = 101
# t = linspace(0, 10, nt)
# p = array(zeros((nt, 3)))
# v = matrix(zeros((9, 1)))
# v[0] = 1
# print v
# for i in xrange(nt):
#     vt, e, hump = spexpv(t[i], out, v)
#     p[i, 0:3] = vt[0:3]
# for i in xrange(3):
#     pl.plot(t, p[:, i], label="(%d)"%(i+1))
# pl.legend(loc='best')
# pl.show()

# dlist = linspace(-5, 5, 201)
# fluo = array([])
# v = matrix(zeros((9, 1)))
# v[0] = 1
# for d in dlist:
#     out = bloch_atomic.copy()
#     out.shift(om13real, laser1_pow)
#     out.shift(d, laser1_det)
#     vt, e, hump = spexpv(200, out, v)
#     fluo = append(fluo, vt[2])
# pl.plot(dlist, fluo)
# pl.show()

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



# ### Time Dependent
out = bloch_atomic.copy()
out.shift(om13real, laser1_pow)
out.shift(d13real, laser1_det)
out.shift(om23real, laser2_pow)
out.shift(d23real, laser2_det)
out.shift(d13real-d23real, laser12_det)

print out


# dd = [0, 4, 8]
# tr2 = {0: 0, 1: 3, 2: 4, 4: 1, 5: 5, 8: 2}
# tr2i = {1:6, 2:7, 5:8}


# tr1 = {0:0, 1:4, 2:8, 3:1, 4:2, 5:5}
# tr1i = {6:1, 7:2, 8:5}
# # conjugates
# tr2c = {1:3, 2:6, 5:7}

# cpm = zeros((9, 9), dtype=complex)

# # Real parts effect on real parts
# for i in tr1.keys():
#     for k in tr1.keys():
#         scale = 1 if (k < 3) else 1/2
#         cpm[tr1[i], tr1[k]] += out[i, k] * scale
# for i in tr1.keys():
#     for k in tr1i.keys():
#         cpm[tr1[i], tr1i[k]] += -1j*out[i, k]/2
#         cpm[tr1[i], tr2c[tr1i[k]]] += 1j*out[i, k]/2
# for i in tr1i.keys():
#     for k in tr1.keys():
#         cpm[tr1i[i], tr1[k]] += 1j*out[i, k]/2
# for i in tr1i.keys():
#     for k in tr1i.keys():
#         cpm[tr1i[i], tr1i[k]] += out[i, k]/2

# for i in tr2.keys():
#     for k in tr2.keys():
#         cpm[i, k] += out[tr2[i], tr2[k]]
#         if tr2i.has_key(k):
#             ipart = 1j*out[tr2[i], tr2i[k]]
#             if i in dd:
#                 ipart /= 2
#             cpm[i, k] += ipart
#             cpm[i, tr2c[k]] += ipart
# for i in tr2.keys():
#     for k in tr2i.keys():
#         print tr2[i], tr2i[k], " -> ", i, k, " = ", out[tr2[i], tr2i[k]]
#         cpm[i, k] += 1j*out[tr2[i], tr2i[k]]

# for i in tr2i.keys():
#     for k in tr2i.keys():
#         cpm[i, k] += 1j*out[tr2i[i], tr2i[k]]
# for i in tr2i.keys():
#     for k in tr2.keys():
#         cpm[i, k] += 1j*out[tr2i[i], tr2[k]]

# Conjugates
# for i in tr2c.keys():
#     cpm[tr2c[i], tr2c[i]] = cpm[i, i].conjugate()

# print real(cpm), "Real"
# print
# print imag(cpm), "Imaginary"
# print out

# for i in xrange(9):
#     v = zeros(9)
#     v[i] = 1
#     w = empty(9)
#     print "%d ----------- " %i
#     print v
#     out.matvec(v, w)
#     print w
#     print


# om13 = om13real
# om23 = om23real
# d12 = d13real - d23real
# d13 = d13real
# d23 = d23real
# obm = array([[0, 0, G31, 0, 0, 0, 0, om13, 0],
#              [0, 0, G32, 0, 0, 0, 0, 0, om23],
#              [0, 0, -Gtot, 0, 0, 0, 0, -om13, -om23],
#              [0, 0, 0, 0, 0, 0, d12, om23/2, om13/2],
#              [0, 0, 0, 0, -Gtot/2, 0, om23/2, d13, 0],
#              [0, 0, 0, 0, 0, -Gtot/2, -om13/2, 0, d23],
#              [0, 0, 0, -d12, -om23/2, om13/2, -Gtot/2, 0, 0],
#              [-om13/2, 0, om13/2, -om23/2, -d13, 0, 0, -Gtot/2, 0],
#              [0, -om23/2, om23/2, -om13/2, 0, -d23, 0, 0, -Gtot/2]])

nt = 151
t = linspace(0, 20,  nt)
p = array(zeros((nt, 9)))
p2 = array(zeros((nt, 9)))
v = matrix(zeros((9, 1)))
v[start] = 1
# print v
for i in xrange(nt):
    vt, e, hump = spexpv(t[i], out, v)
    p[i, 0:9] = vt[0:9]
    # p[i, 0:3] = vt[3:6]
    # v2 = dot(expm(obm*t[i]), v)
    # p2[i, 0:9] = v2[0:9].T
    # p = p2
for i in xrange(3):
    pl.plot(t, p[:, i], '-', label="(%d)"%(i+1))
# for i in xrange(3):
#     pl.plot(t, p2[:, i], 'x', label="X(%d)"%(i+1))
data = loadtxt('../matlab/compare.csv', delimiter=',')
for i in xrange(3):
    pl.plot(t, data[:, i], 'o', label="M(%d)"%(i+1))
pl.legend(loc='best')
pl.title('Populations')

pl.figure()
for i in xrange(3, 6):
    pl.plot(t, p[:, i], label="(%d)"%(i+1))
    pl.plot(t, data[:, i], 'o', label="M(%d)"%(i+1))
pl.title('Real coherence')
pl.legend(loc='best')

pl.figure()
for i in xrange(6, 9):
    pl.plot(t, p[:, i], label="(%d)"%(i+1))
    pl.plot(t, data[:, i], 'o', label="M(%d)"%(i+1))
pl.title('Imaginary coherence')
pl.legend(loc='best')

pl.show()

######################


# print out
# om13 = om13real
# om23 = om23real
# d12 = 0
# d13 = d13real
# d23 = d23real
# obm = array([[0, 0, G31, 0, 0, 0, 0, om13, 0],
#              [0, 0, G32, 0, 0, 0, 0, 0, om23],
#              [0, 0, -Gtot, 0, 0, 0, 0, -om13, -om23],
#              [0, 0, 0, 0, 0, 0, d12, om23/2, om13/2],
#              [0, 0, 0, 0, -Gtot/2, 0, om23/2, d13, 0],
#              [0, 0, 0, 0, 0, -Gtot/2, -om13/2, 0, d23],
#              [0, 0, 0, -d12, -om23/2, om13/2, -Gtot/2, 0, 0],
#              [-om13/2, 0, om13/2, -om23/2, -d13, 0, 0, -Gtot/2, 0],
#              [0, -om23/2, om23/2, -om13/2, 0, -d23, 0, 0, -Gtot/2]])
# print obm

# dlist = linspace(-10, 10, 401)
# fluo = array([])
# v = matrix(zeros((9, 1)))
# v[0] = 1
# for d in dlist:
#     out = bloch_atomic.copy()
#     out.shift(om13real, laser1_pow)
#     out.shift(d, laser1_det)
#     out.shift(om23real, laser2_pow)
#     out.shift(d23real, laser2_det)
#     vt, e, hump = spexpv(100, out, v)
#     fluo = append(fluo, vt[2])
# pl.plot(dlist, fluo, 'r-')
# # data = loadtxt('matlab.csv', delimiter=',')
# # pl.plot(data[:,0], data[:,1], 'k--')
# pl.show()

#     # x11, x22, x33, x12, x13, x23, y12, y13, y23
#     obm = array([[0, 0, G31, 0, 0, 0, 0, om13, 0], \
#                  [0, 0, G32, 0, 0, 0, 0, 0, om23], \
#                  [0, 0, -Gtot, 0, 0, 0, 0, -om13, -om23], \
#                  [0, 0, 0, 0, 0, 0, d12, om23/2, om13/2], \
#                  [0, 0, 0, 0, -Gtot/2, 0, om23/2, d13, 0], \
#                  [0, 0, 0, 0, 0, -Gtot/2, -om13/2, 0, d23], \
#                  [0, 0, 0, -d12, -om23/2, om13/2, -Gtot/2, 0, 0], \
#                  [-om13/2, 0, om13/2, -om23/2, -d13, 0, 0, -Gtot/2, 0], \
#                  [0, -om23/2, om23/2, -om13/2, 0, -d23, 0, 0, -Gtot/2]])
#     for i in range(0, 9):
#         y1[i] += dot(obm[i,:],y)
#     return(y1)





# # Detuning of the two lower levels
# d12 = 0

# # Upper level lifetime, detunings and linewidths
# Gt = 1
# d13 = 0
# G31 = 0.9
# d23 = 0
# G32 = Gt - G31

# # Initial state vector: starting in (1)
# b = array([[1],[0],[0],[0],[0],[0],[0],[0],[0]])

# # Time between pulses
# tpulse = 5
# npulse = 5

# # Total time
# ttotal = tpulse * npulse

# # Strength of pulse on the two transitions
# om0_1 = 2
# om0_2 = 1

# # Pulse "length", square
# om0t = 1

# # Gaussian pulse are calculation
# omgc = om0t/20
# omgb = om0t/2
# omga = om0_1*om0t / (omgc * sqrt(2*pi))

# # put CW area the same as the other two
# om0_1c = om0_1 * om0t / tpulse
# om0_2c = om0_2 * om0t / tpulse

# def omconstant(t):
#     """ Constant omega:
#     warning: global parameters
#     """
#     return (om0_1c, om0_2c)

# def omsquare(t):
#     """ Square pulse omega
#     warning: global parameters
#     """
#     t = mod(t, tpulse)
#     if (t < om0t):
#         return(om0_1, om0_1)
#     else:
#         return(0, 0)

# def gauss(x, a, b, c):
#     """ Gaussian function """
#     return a * exp(-(x-b)**2/(2*c**2))

# def omgauss(t):
#     """ Gaussian pulse
#     warning: global parameters
#     """
#     t = mod(t, tpulse)
#     om = gauss(t, omga, omgb, omgc)
#     return (om, om)

# def f(t0, y, omega):
#     """ Direct integration
#     """
#     y1 = zeros(9)
#     om13, om23 = omega(t0)
#     Gtot = G31 + G32
#     # meaning of d12?
#     d12 = 0
#     # Components of the state vector:
#     # x11, x22, x33, x12, x13, x23, y12, y13, y23
#     # Non-zer0: om31, om32
#     obm = array([[0, 0, G31, 0, 0, 0, 0, om13, 0], \
#                  [0, 0, G32, 0, 0, 0, 0, 0, om23], \
#                  [0, 0, -Gtot, 0, 0, 0, 0, -om13, -om23], \
#                  [0, 0, 0, 0, 0, 0, d12, om23/2, om13/2], \
#                  [0, 0, 0, 0, -Gtot/2, 0, om23/2, d13, 0], \
#                  [0, 0, 0, 0, 0, -Gtot/2, -om13/2, 0, d23], \
#                  [0, 0, 0, -d12, -om23/2, om13/2, -Gtot/2, 0, 0], \
#                  [-om13/2, 0, om13/2, -om23/2, -d13, 0, 0, -Gtot/2, 0], \
#                  [0, -om23/2, om23/2, -om13/2, 0, -d23, 0, 0, -Gtot/2]])
#     for i in range(0, 9):
#         y1[i] += dot(obm[i,:],y)
#     return(y1)


# def doint(y0, t0, dt, tend, omega):
#     """ Do whole length of direct ingegration """
#     r = ode(f).set_integrator('zvode', with_jacobian=False)
#     r.set_initial_value(y0, t0).set_f_params(omega)
#     rt = array([t0])
#     ry = array(y0.T)
#     while r.successful() and r.t < tend:
#         r.integrate(r.t+dt)
#         rt = append(rt, r.t)
#         ry = append(ry, r.y.T, 0)
#     return (rt.T, ry)


# def fmat(t0, dt, omega):
#     """ Bloch matrix exponential at a certain time for a certain
#     interval length.
#     """
#     y1 = zeros(9)
#     om13, om23 = omega(t0+dt/2)
#     Gtot = G31 + G32
#     # meaning of d12?
#     d12 = 0
#     # Components of the state vector:
#     # x11, x22, x33, x12, x13, x23, y12, y13, y23
#     # Non-zer0: om31, om32
#     obm = array([[0, 0, G31, 0, 0, 0, 0, om13, 0], \
#                  [0, 0, G32, 0, 0, 0, 0, 0, om23], \
#                  [0, 0, -Gtot, 0, 0, 0, 0, -om13, -om23], \
#                  [0, 0, 0, 0, 0, 0, d12, om23/2, om13/2], \
#                  [0, 0, 0, 0, -Gtot/2, 0, om23/2, d13, 0], \
#                  [0, 0, 0, 0, 0, -Gtot/2, -om13/2, 0, d23], \
#                  [0, 0, 0, -d12, -om23/2, om13/2, -Gtot/2, 0, 0], \
#                  [-om13/2, 0, om13/2, -om23/2, -d13, 0, 0, -Gtot/2, 0], \
#                  [0, -om23/2, om23/2, -om13/2, 0, -d23, 0, 0, -Gtot/2]])
#     obmt = expm(obm * dt)
# #    for i in range(0, 9):
# #        y1[i] += dot(obmt[i,:],y)
# #    return(y1, obmt)
#     return(obmt)


# def solveexp(y0, tpulse, npulse, omega):
#     """ Do calculation with matrix exponentials for all pulses """
#     obm0 = identity(9)
#     pt = linspace(0, om0t, 41)
#     dt = pt[1] - pt[0]
#     for t in pt[0:-1]:
#         obm0 = dot(fmat(t,dt, omega), obm0)
#     obm0 = dot(fmat(om0t, tpulse-om0t, omega), obm0)
#     y1 = y0.T
#     x = zeros((1, 9))
#     for i in range(0, 9):
#         x[0,i] = dot(obm0[i,:],y1.T)
#     y1 = append(y1, x, 0)
#     for n in range(1, npulse):
#         x = zeros((1, 9))
#         for i in range(0, 9):
#             x[0,i] = dot(obm0[i,:],y1[n,:].T)
#         y1 = append(y1, x, 0)
#     t1 = arange(0, npulse+1) * tpulse
#     return(t1.T, y1)


# # Same pulse area (same power (?)) but CW, Square, Gauss
# dt1 = 0.05
# start = time()
# #rt1, ry1 = doint(b, 0, dt1, ttotal, omconstant)
# #rt1, ry1 = doint(b, 0, dt1, ttotal, omsquare)
# rt1, ry1 = doint(b, 0, dt1, ttotal, omgauss)
# print "Full-numeric: %f s" %(time() - start)


# start = time()
# #te, ye = solveexp(b, tpulse, npulse, omsquare)
# te, ye = solveexp(b, tpulse, npulse, omgauss)
# print "Matrix exponential: %f s" %(time() - start)


# # Plot results
# plot(rt1, ry1[:, 0],'b-', linewidth=2, label='g(1)')
# plot(rt1, ry1[:, 1],'r-', linewidth=2, label='g(2)')
# plot(rt1, ry1[:, 2],'g-', linewidth=2, label='e(3)')
# xlim((0, rt1[-1]))
# for i in range(0, 3):
#     plot(te, ye[:, i], 's')

# legend(loc = 'best')
# show()
