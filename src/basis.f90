subroutine basis(t, nctl, k, u, ind, B)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: basis evaluates the standard b-spline basis
  !               functions at a location x for knot vector t. Adapted
  !               from The NURBS Book, algoritm A2.2

  !     Description of Arguments
  !     Input
  !     t       - Real, Vector of knots. Size (nctl + k)
  !     nctl    - Integer, number of knots
  !     k       - Integer, order of B-spline 
  !     u       - Real, location for evaluation of basis 
  !     ind     - Integer, position in knot vector such that t(ind) <= x
  !
  !     Ouput 
  !     B       - Real, Vector, The k basis vector values
  use precision
  implicit none

  ! Input
  integer,             intent(in)  :: nctl, k, ind
  real(kind=realType), intent(in)  :: t(nctl+k), u

  ! Output
  real(kind=realType), intent(out) :: B(0:k-1)

  ! Working
  real(kind=realType)              :: left(0:k-1), right(0:k-1), temp, saved
  integer                          :: j, r, p

  ! To be consistent with algorithm in The NURBS Book we will use
  ! zero-based ordering here
  B(0) = 1.0
  p = k-1

  do j=1,p
     left(j)  = u - t(ind+1-j)
     right(j) = t(ind+j) - u
     saved = 0.0

     do r=0,j-1
        temp = B(r)/(right(r+1)+left(j-r))
        B(r) = saved+right(r+1)*temp
        saved = left(j-r)*temp
     end do
     
     B(j) = saved
  end do

end subroutine basis

subroutine derivBasis(t, nctl, ku, u, ind, n, Bd)
 
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: derivBasis evaluates the nonzero basis functions and
  !               their derivatives. Algorithm adapted from The NURBS
  !               Book, Algorithm A2.3

  !     Description of Arguments
  !     Input
  !     t       - Real, size(nctl + k) B-spline knot vector
  !     nctl    - Integer, number of control points
  !     ku      - Integer, order of B-spline 
  !     u       - Real, location for basis function evaluations
  !     ind     - Integer, result fro findSpan
  !     n       - Integer, order of derivates to evaluate
  !
  !     Ouput 
  !     Bd      - Real, size(k,k) Basis functions and their derivatives

  use precision
  implicit none
  
  ! Input
  integer,             intent(in)  :: nctl, ku, ind, n
  real(kind=realType), intent(in)  :: t(nctl+ku), u
  
  ! Output
  real(kind=realType), intent(out) :: Bd(0:n,0:n)

  ! Working
  real(kind=realType)              :: left(0:ku-1), right(0:ku-1), saved
  real(kind=realType)              :: a(0:1,0:ku-1), d, temp, ndu(0:ku-1,0:ku-1)
  integer                          :: j, k, r, rk, pk, s1, s2, j1, j2, p

  ! To be consistent with algorithm in The NURBS Book we will use
  ! zero-based ordering here
  ndu(0,0) = 1.0

  ! Define degree as k-1...easier 
  p = ku-1
  do j=1,p
     left(j)  = u - t(ind+1-j)
     right(j) = t(ind+j) - u
     saved = 0.0

     do r=0,j-1
        ndu(j,r) = right(r+1)+left(j-r)
        temp = ndu(r,j-1)/ndu(j,r)

        ndu(r,j) = saved+right(r+1)*temp
        saved = left(j-r)*temp
     end do
        
     ndu(j,j) = saved
  end do

  ! Load the basis functions
  do j=0,p
     Bd(0,j) = ndu(j,p)
  end do

  ! This section computes the derivatives (Eq. [2.9])
  do r=0,p
     s1 = 0  ! Alternate rows in array a
     s2 = 1
     a(0,0) = 1.0
     
     ! Loop to compute kth derivative
     do k=1,n
        d = 0.0
        rk = r-k
        pk = p-k
        
        if (r >= k) then
           a(s2,0) = a(s1,0)/ndu(pk+1,rk)
           d = a(s2,0)*ndu(rk,pk)
        end if
        
        if (rk >= -1) then
           j1 = 1
        else
           j1 = -rk
        end if
        
        if (r-1 <= pk) then
           j2 = k-1
        else
           j2 = p-r
        end if
        
        do j=j1,j2
           a(s2,j) = (a(s1,j)-a(s1,j-1))/ndu(pk+1,rk+j)
           d = d + a(s2,j)*ndu(rk+j,pk)
        end do
        
        if (r <= pk) then
           a(s2,k) = -a(s1,k-1)/ndu(pk+1,r)
           d = d + a(s2,k)*ndu(r,pk)
        end if
        
        Bd(k,r) = d
        
        ! Switch rows
        j = s1
        s1 = s2
        s2 = j
     end do
  end do

  ! Multiply through by the correct factors (Eq 2.9)
  r = p
  do k=1,n
     do j=0,p
        Bd(k,j) = Bd(k,j) * r
     end do
     r = r*(p-k)
  end do

end subroutine derivBasis

subroutine basis_c(t, nctl, k, u, ind, B)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: basis evaluates the standard b-spline basis
  !               functions at a location x for knot vector t. Adapted
  !               from The NURBS Book, algoritm A2.2

  !     Description of Arguments
  !     Input
  !     t       - Real, Vector of knots. Size (nctl + k)
  !     nctl    - Integer, number of knots
  !     k       - Integer, order of B-spline 
  !     u       - Real, location for evaluation of basis 
  !     ind     - Integer, position in knot vector such that t(ind) <= x
  !
  !     Ouput 
  !     B       - Real, Vector, The k basis vector values
  use precision
  implicit none

  ! Input
  integer,             intent(in)  :: nctl, k, ind
  real(kind=realType), intent(in)  :: t(nctl+k)
  complex(kind=realType), intent(in) :: u

  ! Output
  complex(kind=realType), intent(out) :: B(0:k-1)

  ! Working
  complex(kind=realType)              :: left(0:k-1), right(0:k-1), temp, saved
  integer                          :: j, r, p

  ! To be consistent with algorithm in The NURBS Book we will use
  ! zero-based ordering here
  B(0) = cmplx(1.0, 0.0)
  p = k-1

  do j=1,p
     left(j)  = u - t(ind+1-j)
     right(j) = t(ind+j) - u
     saved = 0.0

     do r=0,j-1
        temp = B(r)/(right(r+1)+left(j-r))
        B(r) = saved+right(r+1)*temp
        saved = left(j-r)*temp
     end do
     
     B(j) = saved
  end do

end subroutine basis_c
