   subroutine b2valV(X,Y,N,IDX,IDY,TX,TY,NX,NY,KX,KY,BCOEF,Z)

!***DESCRIPTION
!
!     Written by Gaetan Kenway
!
!     Abstract b2valM is a vectored form of the cmlib function
!         b2val. The intent is provide a vectors X and Y of length N and
!         return a vector Z also of length N . This will provide faster
!         execution when the calling program is Python
!
!     Description of Arguments
!     Input
!     X       - Real, Vector of x-values of size N
!     Y       - Real, Vector of y-values of size N
!     IDX     - Integer,Order of derivative in x. 0 .le. idx .le. kx-1
!                Use 0 for value
!     IDY     - Integer,Order of derivative in y. 0 .le. idy .le. ky-1
!                Use 0 for value
!     TX      - Real,Knot vector for x-direction. Length N+K
!     TY      - Real,Knot vector for y-direction. Length N+K
!     NX      - Integer,Number of control points in x-direction
!     NY      - Integer,Number of control points in y-direction
!     KX      - Integer,order of B-spline in X
!     KY      - Integer,order of B-spline in Y
!     BCOEF   - Real,Array of b-sline coefficients. Size (NX by NY)
!
!     Ouput Z - Real,Array of z-values of size N, interpolated from spline defined
!     by TX,TY,KX,KY, and BCOEF
      
     INTEGER  IDX, IDY, NX, NY, KX, KY, N, I, J
     REAL*8 X(N),Y(N),Z(N)
     REAL*8 BCOEF(NX,NY)      
     REAL*8 TX(NX+KY), TY(NY+KY)
     REAL*8 WORK(3*MAX(KX,KY)+KY)

     external b2val

      do I=1,N
            Z(I) = b2val(X(I),Y(I),IDX,IDY,TX,TY,NX,NY,KX,KY,BCOEF,WORK)
      end do


    end subroutine b2valV
