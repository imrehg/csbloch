from numpy import zeros, float64, array
from expokit import dgpadm

def expf(n, t, A, LDA):
    """to compute:
    exp( t * A(1:n, 1:n) )
    where on input:
    t is a real scalar (positive or negative)
    A(LDA,*) is real n-by-n

    note:
        LDAC is to be provided on input and represent
        the so-called leading-dimension of A.
    """
    ideg = 6;
    lwsp = 4 * n * n + ideg + 1
    wsp = zeros(lwsp, float64)
    iexp = n * n + ideg + 1
    iwsp = zeros(n)
    ns = 0
    iflag = 0
    dgpadm(ideg, t, A, wsp, iwsp, iexp, ns, iflag)
    return wsp[iexp:(iexp+n*n)].reshape((n,n)).T
