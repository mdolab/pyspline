subroutine curve_stiffness(t,n,k,alpha,beta,gpts,K0)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: curve_stiffness generates the stiffness matrix for a curve with a given knot vector t
  !
  !     Description of Arguments
  !     Input
  !     t       - Real, Vector of knots, size n+k
  !     n       - Integer, number of control points
  !     k       - Integer, spline order
  !     alpha   - Real, Stretching stiffness parameter
  !     beta    - Real, Bending stiffness parameter
  !     gpts    - Real, Vector of greville points, size n
  !     Output
  !     K       - Real, Matrix of size nxn, Computed Stiffness

  !    Input
  integer ,           intent(in):: n,k
  double precision,   intent(in):: t(n+k)
  double precision,   intent(in):: alpha,beta
  double precision,   intent(in):: gpts(n)

  !    Output
  double precision,   intent(out):: K0(n,n)
  !    Working
  integer i,j,m,l
  double precision  :: W(3),zeta(3)
  double precision  :: s,factor
  external get_basis

  W(1) = 5.0/9.0
  W(2) = 8.0/9.0
  W(3) = 5.0/9.0

  zeta(1) = -sqrt(0.6)
  zeta(2) = 0
  zeta(3) = sqrt(0.6)

  K0 = 0.0
  do i=1,n
     do j=i,min(i+k,n)
        ! Integrate over the curve
        do m=max(1,i-k),min(n-1,j+k)
           ! Evaluate at Gauss Points
           do l=1,3
              factor = 0.5*(gpts(m+1)-gpts(m))
              s = factor*zeta(l) + 0.5*(gpts(m) + gpts(m+1))
              K0(i,j) = K0(i,j) + w(l)*alpha*get_basis(t,n,k,i,s,1)*get_basis(t,n,k,j,s,1)*factor
              K0(i,j) = K0(i,j) + w(l)*beta*get_basis(t,n,k,i,s,2)*get_basis(t,n,k,j,s,2)*factor
           end do
        end do
     end do
  end do
end subroutine curve_stiffness
