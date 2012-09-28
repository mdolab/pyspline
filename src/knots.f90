subroutine knots_lms(X, N, Nctl, k, t)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract knots_lms generates knots suitable for b-splines
  !
  !     Description of Arguments
  !     Input
  !     X       - Real, size(N) The ordinates to use to make knot vector
  !     N       - Integer, the number of ordinates
  !     Nctl    - Integer, the desired number of control points
  !     k       - order of spline
  !
  !     Ouput 
  !     t       - Real, size(Nctl+k) the resulting knots

  use precision
  implicit none

  ! Input
  integer, intent(in)     :: N
  integer, intent(in)     :: Nctl
  integer, intent(in)     :: K
  real(kind=realType), intent(in)  :: X(N)

  ! Output
  real(kind=realType), intent(out) :: T(Nctl+K)

  ! Working  
  real(kind=realType) d, alpha
  integer  I, J

  ! ----------------------------
  !  PUT K KNOTS AT EACH ENDPOINT -- Knot a knot conditions
  !  ----------------------------
  !
  DO  J=1, K
     t(J) = X(1) ! Left
     t(Nctl+J) = X(N) ! right
  end do

  if (mod(N, 2) == 1) then ! Odd
     d = real(N/(Nctl-k+1.0))
     do j=1, Nctl-k
        i = floor(j*d)
        alpha = j*d-i
        t(k+j) = (1-alpha)*X(i) + alpha*X(i+2)
     end do
  else ! even
     d = real(N/(Nctl-k+1.0))
     do j=1, Nctl-k
        i = floor(j*d)
        alpha = j*d-i+0.5
        t(k+j) = (1-alpha)*X(i) + alpha*X(i+1)
     end do
  end if

end subroutine knots_lms

subroutine knots_interp(X, deriv_ptr, n, nd , k, t)
 !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract knots_lms generates knots suitable for b-spline
  !     interpolation
  !
  !     Description of Arguments
  !     Input
  !     X       - Real, size(N) The ordinates to use to make knot vector
  !     deriv_ptr - Real, size(nd) Flags to determine if a derivative
  !     is specified
  !     n       - Integer, the number of ordinates
  !     nd      - Integer, the number of derivatives specified
  !     k       - order of spline
  !
  !     Ouput 
  !     t       - Real, size(nNctl+k) the resulting knots

  use precision
  implicit none

  ! Input
  integer, intent(in)              :: n, nd, k
  integer, intent(in)              :: deriv_ptr(nd)
  real(kind=realType), intent(in)  :: X(n)

  ! Output
  real(kind=realType), intent(out) :: t(n+nd+k)

  ! Working
  integer                          :: i, j, Nctl

  ! ----------------------------
  !  PUT K KNOTS AT EACH ENDPOINT -- Knot a knot conditions
  !  ----------------------------
  Nctl = n + nd

  DO  J=1, K
     t(J) = X(1) ! Left
     t(Nctl+J) = X(N) ! right
  end do

  if (nd == n) then ! Full Length nd so we can use formulas in NURBS book
     if (k == 3) then
        do i=1, n-2
           t(k+2*i-1) = 0.5*(X(i) + X(i+1))
           t(k+2*i) = X(i+1)
        end do
        t(Nctl) = 0.5*(X(n-1) + 1)
     else if (k ==4) then
        t(5) = X(2)/2
        t(Nctl) = 0.5*(X(n-1) + 1)
         do i=1, n-3
            t(k+2*i  ) = (1.0/3.0)*(2*X(i+1) + X(i+2))
            t(k+2*i+1) = (1.0/3.0)*(X(i+1) + 2*X(i+2))
         end do
     end if

  else if (nd == 0) then
     if (k ==2 ) then
        do j=1, N-k
           t(k+j) = X(j+1)
        end do
     else if (k == 3) then
        do j=1, N-k
           t(k+j) = 0.5*(X(j+1)+X(J+2))
        end do
     else if (k == 4) then
        do j=1, N-k
           t(K+j) = X(j+2)
        end do
     else
        print *, 'Error: Interpolation is only available for k=2, 3 or 4'
        stop
     end if
  else 
     print *, 'Interp_knots with number of derivative != number of points is not yet supported' 
     stop
  end if
end subroutine knots_interp
