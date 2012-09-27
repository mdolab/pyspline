subroutine basis(t, nctl, k, x, ileft, B)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: basis evaluates the standard b-spline basis functions at 
  !               a location x for knot vector t

  !     Description of Arguments
  !     Input
  !     t       - Real, Vector of knots. Size (nctl + k)
  !     nctl    - Integer, number of knots
  !     k       - Integer, order of B-spline 
  !     x       - Real, location for evaluation of basis 
  !     ileft   - Integer, position in knot vector such that t(ileft) <= x
  !
  !     Ouput 
  !     B       - Real, Vector, The k basis vector values
  use precision
  implicit none

  ! Input
  real(kind=realType), intent(in)  :: t(nctl+k)
  integer,             intent(in)  :: nctl, k, ileft

  ! Output
  real(kind=realType), intent(out) :: B(k)

  ! Working
  real(kind=realType)              :: vm, vmprev, x, diff1(k), diff2(k)
  integer                          ::  l, j, ipj, jp1ml, jp1
    
  B(:) = 0.0
  B(1) = 1.0

  do j=1, k-1
     ipj = ileft + j
     diff1(j) = t(ipj)-x
     diff2(j) = x-t(ileft-j+1)
     vmprev = 0.0
     jp1 = j+1
     do l=1, j
        jp1ml = jp1-l
        vm = B(l)/(diff1(l)+diff2(jp1ml))
        B(l) = vm*diff1(l) + vmprev
        vmprev = vm*diff2(jp1ml)
     end do
     B(jp1) = vmprev
  end do

end subroutine basis

subroutine derivBasis(ind, u, ku, n, nctl, t, Bd)
 
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: derivBasis evaluates the nonzero basis functions and
  !               their derivatives. Algorithm adapted from The NURBS
  !               Book, Algorithm A2.3

  !     Description of Arguments
  !     Input
  !     ind     - Integer, result fro findSpan
  !     u       - Real, location for basis function evaluations
  !     k       - Integer, order of B-spline 
  !     n       - Integer, order of derivates to evaluate
  !     nctl    - Integer, number of control points
  !     t       - Real, size(nctl + k) B-spline knot vector
  !
  !     Ouput 
  !     Bd      - Real, size(k,k) Basis functions and their derivatives

  use precision
  implicit none
  
  ! Input
  real(kind=realType), intent(in)  :: t(nctl+ku), u
  integer,             intent(in)  :: nctl, ku, ind, n
  
  ! Output
  real(kind=realType), intent(out) :: Bd(0:ku-1,0:ku-1)

  ! Working
  real(kind=realType)              :: ndu(0:ku-1,0:ku-1), left(0:ku-1), right(0:ku-1), saved
  real(kind=realType)              :: a(0:1,0:ku-1), d, temp
  integer                          :: j, k, r, rk, pk, s1, s2, j1, j2, p

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

  ! This section computes the derivates (Eq. [2.9])

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
