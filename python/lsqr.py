"""
Solve the least-squares problem

  minimize ||Ax-b||

using LSQR.  This is a line-by-line translation from Matlab code
available at http://www.stanford.edu/~saunders/lsqr.

Michael P. Friedlander, University of British Columbia
Dominique Orban, Ecole Polytechnique de Montreal
$Id: lsqr.py 88 2008-09-29 04:43:15Z d-orban $
"""

from numpy import zeros, dot
from numpy.linalg import norm
from math import sqrt

# Simple shortcuts---linalg.norm is too slow for small vectors
def normof2(x,y): return sqrt(x**2 + y**2)
def normof4(x1,x2,x3,x4): return sqrt(x1**2 + x2**2 + x3**2 + x4**2)

def lsqr( m, n, aprod, b, damp, atol, btol, conlim, itnlim, show, wantvar = False ):
    """
    [ x, istop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm, var ]...
    = lsqr( m, n, @aprod, b, damp, atol, btol, conlim, itnlim, show );
    
    LSQR solves  Ax = b  or  min ||b - Ax||_2  if damp = 0,
              or   min || (b)  -  (  A   )x ||   otherwise.
                       || (0)     (damp I)  ||2
    A  is an m by n matrix defined by  y = aprod( mode,m,n,x ),
    where @Aprod refers to a function handle that
    performs the matrix-vector operations.
    If mode = 1,   aprod  must return  y = Ax   without altering x.
    If mode = 2,   aprod  must return  y = A'x  without altering x.
    WARNING:   The file containing the function 'aprod'
               must not be called aprodname.m !!!!

    ----------------------------------------------------------------------
    LSQR uses an iterative (conjugate-gradient-like) method.
    For further information, see 
    1. C. C. Paige and M. A. Saunders (1982a).
       LSQR: An algorithm for sparse linear equations and sparse least squares,
       ACM TOMS 8(1), 43-71.
    2. C. C. Paige and M. A. Saunders (1982b).
       Algorithm 583.  LSQR: Sparse linear equations and least squares problems,
       ACM TOMS 8(2), 195-209.
    3. M. A. Saunders (1995).  Solution of sparse rectangular systems using
       LSQR and CRAIG, BIT 35, 588-604.

    Input parameters:
    iw, rw      are not used by lsqr, but are passed to aprod.
    atol, btol  are stopping tolerances.  If both are 1.0e-9 (say),
                the final residual norm should be accurate to about 9 digits.
                (The final x will usually have fewer correct digits,
                depending on cond(A) and the size of damp.)
    conlim      is also a stopping tolerance.  lsqr terminates if an estimate
                of cond(A) exceeds conlim.  For compatible systems Ax = b,
                conlim could be as large as 1.0e+12 (say).  For least-squares
                problems, conlim should be less than 1.0e+8.
                Maximum precision can be obtained by setting
                atol = btol = conlim = zero, but the number of iterations
                may then be excessive.
    itnlim      is an explicit limit on iterations (for safety).
                show = 1    gives an iteration log,
                show = 0    suppresses output.

    Output parameters:
    x           is the final solution.
    istop       gives the reason for termination.
    istop       = 1 means x is an approximate solution to Ax = b.
                = 2 means x approximately solves the least-squares problem.
                r1norm      = norm(r), where r = b - Ax.
                r2norm      = sqrt( norm(r)^2  +  damp^2 * norm(x)^2 )
                = r1norm if damp = 0.
    anorm       = estimate of Frobenius norm of Abar = [  A   ].
                                                       [damp*I]
    acond       = estimate of cond(Abar).
    arnorm      = estimate of norm(A'*r - damp^2*x).
    xnorm       = norm(x).
    var         (if present) estimates all diagonals of (A'A)^{-1} (if damp=0)
                or more generally (A'A + damp^2*I)^{-1}.
                This is well defined if A has full column rank or damp > 0.
                (Not sure what var means if rank(A) < n and damp = 0.)
                
    1990: Derived from Fortran 77 version of LSQR.
    22 May 1992: bbnorm was used incorrectly.  Replaced by anorm.
    26 Oct 1992: More input and output parameters added.
    01 Sep 1994: Matrix-vector routine is now a parameter 'aprodname'.
                 Print log reformatted.
    14 Jun 1997: show  added to allow printing or not.
    30 Jun 1997: var   added as an optional output parameter.
    07 Aug 2002: Output parameter rnorm replaced by r1norm and r2norm.
                 Michael Saunders, Systems Optimization Laboratory,
                 Dept of MS&E, Stanford University.
    08 Dec 2007: Line-by-line translation to Python
    ----------------------------------------------------------------------
    """

    # Initialize.

    msg=['The exact solution is  x = 0                              ',
         'Ax - b is small enough, given atol, btol                  ',
         'The least-squares solution is good enough, given atol     ',
         'The estimate of cond(Abar) has exceeded conlim            ',
         'Ax - b is small enough for this machine                   ',
         'The least-squares solution is good enough for this machine',
         'Cond(Abar) seems to be too large for this machine         ',
         'The iteration limit has been reached                      ']

    if wantvar:
        var = zeros(n,1)
    else:
        var = None
    
    if show:
        print ' '
        print 'LSQR            Least-squares solution of  Ax = b'
        str1 = 'The matrix A has %8g rows  and %8g cols' % (m, n)
        str2 = 'damp = %20.14e     wantvar = %-5s' % (damp, repr(wantvar))
        str3 = 'atol = %8.2e                 conlim = %8.2e' % (atol, conlim)
        str4 = 'btol = %8.2e                 itnlim = %8g' % (btol, itnlim)
        print str1;   print str2;   print str3;   print str4;
        
    itn    = 0;   istop  = 0;   nstop  = 0
    ctol   = 0.0
    if conlim > 0.0: ctol = 1.0/conlim
    anorm  = 0.;	acond  = 0.
    dampsq = damp**2;	ddnorm = 0.;		res2   = 0.
    xnorm  = 0.;	xxnorm = 0.;		z      = 0.
    cs2    = -1.;	sn2    = 0.
    
    # Set up the first vectors u and v for the bidiagonalization.
    # These satisfy  beta*u = b,  alfa*v = A'u.
    
    u      = b[:m];	x    = zeros(n)
    alfa   = 0.;	beta = norm( u )
    if beta > 0:
        u = (1.0/beta) * u;	v = aprod(2, m, n, u)
        alfa = norm( v );

    if alfa > 0:
        v = (1.0/alfa) * v;    w = v.copy();
        
    arnorm = alfa * beta;
    if arnorm == 0:
        print msg[0]
        return (x, istop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm, var)

    rhobar = alfa;		phibar = beta;		bnorm  = beta;
    rnorm  = beta
    r1norm = rnorm
    r2norm = rnorm
    head1  = '   Itn      x(1)       r1norm     r2norm '
    head2  = ' Compatible   LS      Norm A   Cond A'
    
    if show:
        print ' '
        print head1+head2
        test1  = 1.0;		test2  = alfa / beta
        str1   = '%6g %12.5e'     % (itn,   x[0])
        str2   = ' %10.3e %10.3e' % (r1norm, r2norm)
        str3   = '  %8.1e %8.1e'  % (test1,  test2)
        print str1+str2+str3
        
    # ------------------------------------------------------------------
    #     Main iteration loop.
    # ------------------------------------------------------------------
    while itn < itnlim:
        itn = itn + 1
        #   Perform the next step of the bidiagonalization to obtain the
        #   next  beta, u, alfa, v.  These satisfy the relations
        #              beta*u  =  a*v   -  alfa*u,
        #              alfa*v  =  A'*u  -  beta*v.

        u    = aprod(1, m, n, v)  -  alfa*u
        beta = norm( u );
        if beta > 0:
            u     = (1.0/beta) * u
            anorm = normof4(anorm, alfa, beta, damp)
            v     = aprod(2, m, n, u)  -  beta*v
            alfa  = norm( v )
            if alfa > 0:  v = (1.0/alfa) * v
        
        # Use a plane rotation to eliminate the damping parameter.
        # This alters the diagonal (rhobar) of the lower-bidiagonal matrix.

        rhobar1 = normof2(rhobar, damp)
        cs1     = rhobar / rhobar1
        sn1     = damp   / rhobar1
        psi     = sn1 * phibar
        phibar  = cs1 * phibar
        
        #  Use a plane rotation to eliminate the subdiagonal element (beta)
        # of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
        
        rho     =   normof2(rhobar1, beta)
        cs      =   rhobar1/ rho
        sn      =   beta   / rho
        theta   =   sn * alfa
        rhobar  = - cs * alfa
        phi     =   cs * phibar
        phibar  =   sn * phibar
        tau     =   sn * phi
        
        # Update x and w.
            
        t1      =   phi  /rho;
        t2      = - theta/rho;
        dk      =   (1.0/rho)*w;
            
        x       = x      +  t1*w
        w       = v      +  t2*w
        ddnorm  = ddnorm +  norm(dk)**2
        if wantvar: var = var  +  dk*dk
            
        # Use a plane rotation on the right to eliminate the
        # super-diagonal element (theta) of the upper-bidiagonal matrix.
        # Then use the result to estimate  norm(x).
            
        delta   =   sn2 * rho
        gambar  = - cs2 * rho
        rhs     =   phi  -  delta * z
        zbar    =   rhs / gambar
        xnorm   =   sqrt(xxnorm + zbar**2)
        gamma   =   normof2(gambar, theta)
        cs2     =   gambar / gamma
        sn2     =   theta  / gamma
        z       =   rhs    / gamma
        xxnorm  =   xxnorm  +  z**2
            
        # Test for convergence.
        # First, estimate the condition of the matrix  Abar,
        # and the norms of  rbar  and  Abar'rbar.
        
        acond   =   anorm * sqrt( ddnorm )
        res1    =   phibar**2
        res2    =   res2  +  psi**2
        rnorm   =   sqrt( res1 + res2 )
        arnorm  =   alfa * abs( tau )
        
        # 07 Aug 2002:
        # Distinguish between
        #    r1norm = ||b - Ax|| and
        #    r2norm = rnorm in current code
        #           = sqrt(r1norm^2 + damp^2*||x||^2).
        #    Estimate r1norm from
        #    r1norm = sqrt(r2norm^2 - damp^2*||x||^2).
        # Although there is cancellation, it might be accurate enough.
        
        r1sq    =   rnorm**2  -  dampsq * xxnorm
        r1norm  =   sqrt( abs(r1sq) )
        if r1sq < 0: r1norm = - r1norm
        r2norm  =   rnorm
        
        # Now use these norms to estimate certain other quantities,
        # some of which will be small near a solution.
        
        test1   =   rnorm / bnorm
        test2   =   arnorm/( anorm * rnorm )
        test3   =       1.0 / acond
        t1      =   test1 / (1    +  anorm * xnorm / bnorm)
        rtol    =   btol  +  atol *  anorm * xnorm / bnorm
        
        # The following tests guard against extremely small values of
        # atol, btol  or  ctol.  (The user may have set any or all of
        # the parameters  atol, btol, conlim  to 0.)
        # The effect is equivalent to the normal tests using
        # atol = eps,  btol = eps,  conlim = 1/eps.
        
        if itn >= itnlim:   istop = 7
        if 1 + test3  <= 1: istop = 6
        if 1 + test2  <= 1: istop = 5
        if 1 + t1     <= 1: istop = 4
        
        # Allow for tolerances set by the user.
        
        if  test3 <= ctol:  istop = 3
        if  test2 <= atol:  istop = 2
        if  test1 <= rtol:  istop = 1
        
        # See if it is time to print something.
            
        prnt = False;
        if n     <= 40       : prnt = True
        if itn   <= 10       : prnt = True
        if itn   >= itnlim-10: prnt = True
        if itn % 10 == 0     : prnt = True
        if test3 <=  2*ctol  : prnt = True
        if test2 <= 10*atol  : prnt = True
        if test1 <= 10*rtol  : prnt = True
        if istop !=  0       : prnt = True
        
        if prnt:
            if show:
                str1 = '%6g %12.5e'     %(   itn,   x[0] )
                str2 = ' %10.3e %10.3e' %(r1norm, r2norm )
                str3 = '  %8.1e %8.1e'  %( test1,  test2 )
                str4 = ' %8.1e %8.1e'   %( anorm,  acond )
                print str1+str2+str3+str4
                
        if istop > 0: break

    # End of iteration loop.
    # Print the stopping condition.
        
    if show:
        print ' '
        print 'LSQR finished'
        print msg[istop]
        print ' '
        str1 = 'istop =%8g   r1norm =%8.1e'   %(istop, r1norm )
        str2 = 'anorm =%8.1e   arnorm =%8.1e' %(anorm, arnorm )
        str3 = 'itn   =%8g   r2norm =%8.1e'   %(  itn, r2norm )
        str4 = 'acond =%8.1e   xnorm  =%8.1e' %(acond, xnorm  )
        str5 = '                  bnorm  =%8.1e'    % bnorm
        print str1 + '   ' + str2
        print str3 + '   ' + str4
        print str5
        print ' '
        
    return ( x, istop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm, var )
