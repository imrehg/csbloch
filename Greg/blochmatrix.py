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
