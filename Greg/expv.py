from __future__ import division
from numpy import *
from scipy.linalg import expm
from pysparse import spmatrix

def norminf(A):
    """ Infinite norm, equivalent to Matlab's norm(A, 'inf') """
    return float(max(abs(A.sum(axis=1))))

def norm(v):
    """ L2 norm, equivalent to Matlab's norm(v) """
    return float(sqrt(dot(v[:].T,v[:])))

def expv(t, A, v, m = 30):
    A = matrix(A)
    n = A.shape[0]
    tol = 1e-7
    m = min([n, m])
    anorm = norminf(A)
    
    mxrej = 10
    btol  = 1.0e-7
    gamma = 0.9
    delta = 1.2; 
    mb    = m
    t_out   = abs(t);
    nstep = 0
    t_new   = 0
    t_now = 0
    s_error = 0;
    rndoff= anorm * finfo(float).eps

    k1 = 2
    xm = 1/m
    normv = norm(v)
    beta = normv
    fact = (((m+1)/exp(1))**(m+1))*sqrt(2*pi*(m+1))
    t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))**xm
    s = 10**(floor(log10(t_new))-1)
    t_new = ceil(t_new/s)*s
    sgn = sign(t)
    nstep = 0

    w = v.copy();
    hump = normv;
    while t_now < t_out:
        nstep += 1
        t_step = min( t_out-t_now, t_new )
        V = zeros((n, m+2), float64)
        H = zeros((m+2, m+2), float64)

        V[0:n, 0] = (1/beta)*w.T.copy()
        for j in xrange(m):
            p = dot(A, V[:, j])
            # p =  (A * V[:, j])
            for i in xrange(j+1):
                H[i, j] = dot(V[:,i], p.T)
                p = p - H[i,j] * V[:,i]
            s = norm(p.T)
            if s < btol:
                   k1 = 0
                   mb = j
                   t_step = t_out-t_now
                   break
            H[j+1,j] = s;
            V[0:n, j+1] = (1/s)*p.copy()

        if k1 != 0: 
            H[m+1,m] = 1;
            avnorm = norm(dot(A,V[:, m]).T)
        ireject = 0
        while ireject <= mxrej:
            mx = mb + k1 + 1
            F = expm(sgn * t_step * H[:mx, :mx])
            if k1 == 0:
                err_loc = btol
                break
            else:
                phi1 = abs( beta*F[m,1] );
                phi2 = abs( beta*F[m+1,1] * avnorm );
                if phi1 > 10*phi2:
                    err_loc = phi2
                    xm = 1/m
                elif phi1 > phi2:
                    err_loc = (phi1*phi2)/(phi1-phi2)
                    xm = 1/m;
                else:
                    err_loc = phi1
                    xm = 1/(m-1)
                
            if err_loc <= delta * t_step*tol:
                break
            else:
                t_step = gamma * t_step * (t_step*tol/err_loc)**xm
                s = 10**(floor(log10(t_step))-1)
                t_step = ceil(t_step/s) * s
                if ireject == mxrej:
                    print 'The requested tolerance is too high.'
                ireject += 1

        mx = mb + max(0, k1-1) + 1
        w = dot(V[:, :mx], beta*F[:mx, 0])
        beta = norm(w)
        hump = max(hump, beta)

        t_now = t_now + t_step;
        t_new = gamma * t_step * (t_step*tol/err_loc)**xm;
        s = 10**(floor(log10(t_new))-1)
        t_new = ceil(t_new/s) * s

        err_loc = max(err_loc,rndoff)
        s_error = s_error + err_loc

    err = s_error
    hump = hump / normv
    return matrix(w).T, err, hump

def spexpv(t, A, v, tol=1e-7, m = 30):
    """
    Sparse matrix exponentiation using the Krylov subspace method.
    w, err, hump = spexpv(t, A, v, tol=1e-7, m = 30)

    Not calculating the matrix exponential directly but it's effect on a vector:
    w = expm(tA)v

    Input:
    ======
    t : time (can be t<0)
    A : propagation matrix (in pysparse.spmatrix.ll_mat format)
    v : initial state vector (in numpy.array format)
    tol : tolerance
    m : maximum number of used eigenvalues (if A has dimension larger than m)

    Output:
    =======
    w : exp(tA)v vector, in the same format as scipy.linalg.expm
    err : final error (should be err < tol, or very close to that)
    hump : max(||exp(sA)||) for s in [0, t], or s in [t, 0] if  t < 0
           It measures the conditioning of the problem. It is good if hump = 1,
           wheres poor if hump >> 1.

    The Krylov-subspace calculation has a number of different possible implementations,
    here the Arnoldi-iteration is used.

    Code is based on expokit and in particular the Matlab implementation,
    http://www.maths.uq.edu.au/expokit/matlab/expv.m
    Roger B. Sidje (rbs@maths.uq.edu.au)
    EXPOKIT: Software Package for Computing Matrix Exponentials.
    ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
    """
    n = A.shape[0]
    m = min([n, m])
    anorm = A.norm('inf')

    mxrej = 10
    btol  = 1.0e-7
    gamma = 0.9
    delta = 1.2;
    mb    = m
    t_out   = abs(t);
    nstep = 0
    t_new   = 0
    t_now = 0
    s_error = 0;
    rndoff= anorm * finfo(float).eps

    k1 = 2
    xm = 1/m
    normv = norm(v)
    beta = normv
    fact = (((m+1)/exp(1))**(m+1))*sqrt(2*pi*(m+1))
    t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))**xm
    s = 10**(floor(log10(t_new))-1)
    t_new = ceil(t_new/s)*s
    sgn = sign(t)
    nstep = 0

    w = v.copy();
    hump = normv;
    # Arnoldi iteration
    while t_now < t_out:
        nstep += 1
        t_step = min( t_out-t_now, t_new )
        V = zeros((n, m+2), float64)
        H = zeros((m+2, m+2), float64)

        V[0:n, 0] = (1/beta)*w.T.copy()
        for j in xrange(m):
            p = array(zeros(n), dtype=float64)
            A.matvec(V[:, j], p)
            for i in xrange(j+1):
                H[i, j] = dot(V[:,i], p.T)
                p = p - H[i,j] * V[:,i]
            s = norm(p.T)
            if s < btol:
                   k1 = 0
                   mb = j
                   t_step = t_out-t_now
                   break
            H[j+1,j] = s;
            V[0:n, j+1] = (1/s)*p.copy()

        # If tolerance level is above threshold (mostly if m < n)
        if k1 != 0:
            H[m+1,m] = 1;
            res = array(zeros(n), dtype=float64)
            A.matvec(V[:, m], res)
            avnorm = norm(res)
        ireject = 0
        while ireject <= mxrej:
            mx = mb + k1 + 1
            F = expm(sgn * t_step * H[:mx, :mx])
            if k1 == 0:
                err_loc = btol
                break
            else:
                phi1 = abs( beta*F[m,1] );
                phi2 = abs( beta*F[m+1,1] * avnorm );
                if phi1 > 10*phi2:
                    err_loc = phi2
                    xm = 1/m
                elif phi1 > phi2:
                    err_loc = (phi1*phi2)/(phi1-phi2)
                    xm = 1/m;
                else:
                    err_loc = phi1
                    xm = 1/(m-1)

            if err_loc <= delta * t_step*tol:
                break
            else:
                t_step = gamma * t_step * (t_step*tol/err_loc)**xm
                s = 10**(floor(log10(t_step))-1)
                t_step = ceil(t_step/s) * s
                if ireject == mxrej:
                    print 'The requested tolerance is too high.'
                ireject += 1

        mx = mb + max(0, k1-1) + 1
        w = dot(V[:, :mx], beta*F[:mx, 0])
        beta = norm(w)
        hump = max(hump, beta)

        t_now = t_now + t_step;
        t_new = gamma * t_step * (t_step*tol/err_loc)**xm;
        s = 10**(floor(log10(t_new))-1)
        t_new = ceil(t_new/s) * s

        err_loc = max(err_loc,rndoff)
        s_error = s_error + err_loc

    err = s_error
    hump = hump / normv
    return matrix(w).T, err, hump

if __name__ == "__main__":
    from time import time
    rep = 5
    n = 150
    A = random.rand(n, n)
    # A = matrix([[1., 5., 3.], [-1., 0., 1.], [3., 2., 1.]])

    for i in xrange(n):
        for j in xrange(n):
            if (random.randn() > 0.001):
                A[i, j] = 0.

    v = eye(n, 1)
    t = 1
    start = time()
    for r in range(rep):
        w, err, hump = expv(t,A,v)
    print "Krylov       :", time()-start

    As = spmatrix.ll_mat(n, n)
    j = range(n)
    for i in xrange(n):
        for j in xrange(n):
            As[i, j] = A[i, j]

    start = time()
    for r in range(rep):
        w, err, hump = spexpv(t, As, v)
    print "Sparse Krylov:", time() - start

    start = time()
    for i in xrange(rep):
        w2 = dot(expm(t*A),v)
    print "Scipy        :", time() - start
    # print (w2 - w)/w
