function get_basis(t,n,k,i,s,ideriv)

!***DESCRIPTION
!
!     Written by Gaetan Kenway
!
!     Abstract: Get the b-spline basis function value for
!     knots, t, and interval, i.  
!
!     Description of Arguments
!     Input
!     t       - Real, Vector of knots, size n+k
!     n       - Integer, number of control points
!     k       - Integer, spline order
!     i       - Integer, the i-th basis function
!     s       - Real, the location for spline basis eval
!     ideriv  - Integer, the derivative order, 0 for val, 1 for first derivative ect.
!
!     Output
!     val      - Real, the value of the request basis function
 
!    Input
     INTEGER n,k,i,ideriv
     REAL*8 t(n+k),s
!    Output
     REAL*8 get_basis
!    Working
     real*8 coef(n),work(3*k)
     integer inbv
     
     external bvalu
     inbv = 1
     coef(:) = 0.0
     coef(i) = 1.0

     get_basis = bvalu(t,coef,n,k,ideriv,s,inbv,work)
   end function get_basis
