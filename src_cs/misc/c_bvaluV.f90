      subroutine bvaluV(T,A,N,K,IDERIV,X,Y,NX,INBV,WORK)
!***DESCRIPTION
!
!     Written by Gaetan Kenway
!
!     Abstract
!         bvaluV is a vectored form of the cmlib function bvalu
!         The intent is provide a vector X of length NX and return
!         a vector of evaluated values Y (of length NX). This will
!         provide faster execution when the calling program is Python
!
!     Description of Arguments
!         Input
!          T       - knot vector of length N+K
!          A       - B-spline coefficient vector of length N
!          N       - number of B-spline coefficients
!                    N = sum of knot multiplicities-K
!          K       - order of the B-spline, K .GE. 1
!          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
!                    IDERIV=0 returns the B-spline value
!          X       - input vector of length NX. T(K) .LE. X(I) .LE. T(N+1)
!          Y       - output vector of lenght NX. 
!          NX      - integer. Length of array X
!          INBV    - an initialization parameter which must be set
!                    to 1 the first time BVALU is called.
!
!         Output
!          INBV    - INBV contains information for efficient process-
!                    ing after the initial call and INBV must not
!                    be changed by the user.  Distinct splines require
!                    distinct INBV parameters.
!          WORK    - work vector of length 3*K.
!          BVALU   - value of the IDERIV-th derivative at X
!
      
        use complexify
	implicit complex*16(a-h, o-z) 
      INTEGER IDERIV,N,K,NX,I
      
      complex*16 A(N)
      complex*16 T(N+K)
      complex*16 X(NX)
      complex*16 Y(NX)
      complex*16 WORK(3*K)

      EXTERNAL bvalu

      do I = 1,NX
         if (X(I) .gt. T(N+K)) then
            X(I) = T(N+K)
         end if
         if (X(I) .lt. T(1)) then
            X(I) = T(1)
         end if


         Y(I)= BVALU(T,A,N,K,IDERIV,X(I),INBV,WORK)         
      end do
      
    end subroutine bvaluV
