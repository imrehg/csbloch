from __future__ import division
import pysparse.spmatrix as sp
import numpy as np
import physicspy.quantum as quantum

class HFSLevel:
    """ A hyperfine level """

    def __init__(self, L, J, Ahfs, Bhfs, Chfs, gJ, G, name=""):
        self.params = {'L': L,
                       'J': J,
                       'Ahfs': Ahfs,
                       'Bhfs': Bhfs,
                       'Chfs': Chfs,
                       'gJ': gJ,
                       'G': G,
                       'name': name,
                       }

    def quantum(self):
        return self.params

class Level:
    """ A level with all its parameters """

    L, J, F, mF, hfs = None, None, None, None, None

    def __init__(self, L, J, F, mF, hfs):
        """ Store values """
        self.L = L
        self.J = J
        self.F = F
        self.mF = mF
        pass

    def __repr__(self):
        """ How to print it """
        return("Level: L=%g J=%g F=%g mF=%g\n" %(self.L, self.J, self.F, self.mF))


class Atom:
    """ A collection of levels and reference between them """

    I = 0
    levels = []

    def __init__(self, I, gI, gS, gL, uB, lower, upper):
        self.I = I
        self.gI = gI
        self.gS = gS
        self.gL = gL
        self.uB = uB
        self.lower = lower
        self.upper = upper

        for hfslevel in [self.lower, self.upper]:
            lpar = hfslevel.quantum()
            # Need to turn fractional, that is float here, into integer for range
            F = range(int(np.abs(lpar['J']-self.I)), int(lpar['J']+self.I+1))
            for Fi in F:
                for mF in xrange(-Fi, Fi+1):
                    self.levels += [Level(lpar['L'], lpar['J'], Fi, mF, hfslevel)]
        self.nlevel = len(self.levels)


    def bloch(self):
        """ Return all the bloch matrices """
        return (self.mblochbase(),
                self.mmagfield(),
                self.mdetuning(),
                self.mlaserpi(),
                self.mlasersp(),
                self.mlasersm(),
                )

    ### Can these be combined not to duplicate some calculations?
    def mblochbase(self):
        """ Base bloch matrix """
        bbase = sp.ll_mat(self.nlevel**2, self.nlevel**2)
        pass

    def mmagfield(self):
        """ Bloch matrix of level shifts for unit magnetic field """
        pass

    def mdetuning(self):
        """ Bloch matrix for detunings """
        pass

    def mlaserpi(self):
        """ Interacton terms of bloch matrix for pi-polarized light """
        pass

    def mlasersp(self):
        """ Interacton terms of bloch matrix for sigma-plus-polarized light """
        pass

    def mlasersm(self):
        """ Interacton terms of bloch matrix for sigma-minus-polarized light """
        pass


### Constants and parameters for Cs
I = 7/2
gI = -0.00039885395
gS = 2.0023193043622
gL = 0.99999587
uB = 1.399624e6

J = np.array([1/2, 3/2])
Ahfs = np.array([2.2981572425e9, 50.28827e6])
Bhfs = np.array([0, -0.4934e6])
Chfs = np.array([0, 0.560e3])
gJ = np.array([2.00254032, 1.334])

### The levels we are interested in for the D2 line
s12 = HFSLevel(0, 1/2, Ahfs[0], Bhfs[0], Chfs[0], gJ[0], G=0, name="S1/2")
p32 = HFSLevel(1, 3/2, Ahfs[1], Bhfs[1], Chfs[1], gJ[1], G=1, name="P3/2")

Cs = Atom(7/2, gI, gS, gL, uB, lower=s12, upper=p32)
Cs.bloch()

# def detuning(qconst, qnumber):
#     Ahfs, Bhfs, Chfs, gJ = qconst
#     I, J, F = qnumber
#     energy = quantum.energyhfsalkali(F,J,I,Ahfs,Bhfs,Chfs)

#     gF = gJ * ( F*(F+1) - I*(I+1) + J*(J+1) ) / ( 2*F*(F+1) ) + \
#          gI * ( F*(F+1) + I*(I+1) - J*(J+1) ) / ( 2*F*(F+1) )
#     print gF
#     print("Detuning: MHz/ G : %g" %(uB*gF/1e6))

#     print ((gJ - gI)**2 * uB**2) / 2

#     return energy


# i = 1
# F = 3
# qconst = (Ahfs[i], Bhfs[i], Chfs[i], gJ[i])
# qnumber = (I, J[i], F)
# print detuning(qconst, qnumber)

# # # The full-level calculation
# # # 48**2 x 48**2 matrix
# # bloch_base = ll_mat(2304, 2304, 1000)

# # # No upper level coherence due to collisions (very short lifetime)
# # # (16**2 + 32) x (16**2 + 32)
# # bloch_base = ll_mat(288, 288, 1000)


