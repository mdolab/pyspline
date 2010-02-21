subroutine knots_interp(X,deriv_ptr,n,nd,k,t)
  !  --------------------------------------------------------------------
  !  knots chooses a knot sequence,T, for interpolation
  !  --------------------------------------------------------------------

  implicit none
  integer, intent(in)     :: n,nd,k
  integer, intent(in)     :: deriv_ptr(nd)
  double precision, intent(in)  :: X(n) ! ie s values
  double precision, intent(out) :: T(n+nd+k)
  !
  ! LOCAL VARIABLES
  !
  double precision d,alpha,temp(n-2)
  INTEGER  I, J, Nctl,d_count
  ! ----------------------------
  !  PUT K KNOTS AT EACH ENDPOINT -- Knot a knot conditions
  !  ----------------------------
  Nctl = n + nd

  DO  J=1,K
     T(J) = X(1) ! Left
     T(Nctl+J) = X(N) ! right
  end do

  if (nd == n) then ! Full Length nd so we can use formulas in NURBS book
     if (k == 3) then
        do i=1,n-2
           T(k+2*i-1) = 0.5*(X(i) + X(i+1))
           T(k+2*i) = X(i+1)
        end do
        T(Nctl) = 0.5*(X(n-1) + 1)
     else if (k ==4) then
        T(5) = X(2)/2
        T(Nctl) = 0.5*(X(n-1) + 1)
         do i=1,n-3
            T(k+2*i  ) = (1.0/3.0)*(2*X(i+1) + X(i+2))
            T(k+2*i+1) = (1.0/3.0)*(X(i+1) + 2*X(i+2))
         end do
     end if

  else if (nd == 0) then
     do j=1,N-k
        T(K+j) = X(j+2)
     end do
  else
     print *,'Interp_knots with number of derivative != number of points is not yet supported'
     stop
  end if
end subroutine knots_interp
