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

  ! Do a "regular" interpolation knot vector in 'temp' and then add in
  ! then additional ones we need for the derivative constraints
  !do j = 1,N-2 ! Standard cubic spline interpolation (points at knots)
  !   temp(j) = X(j+1)
  !end do
  !print *,'X:',X
  !print *,'N:',N
  !print *,'temp:',temp
  ! if (nd > 0) then ! we have to add additional knots 
  !    print *,'in here'
  !    print *,'deriv_ptr:',deriv_ptr
  !    d_count = -1
  !    do j=2,N-1
        
  !       if (deriv_ptr(j) +1 == j) then !
  !          print *,'j:',j
  !          print *,'setting indices:',k+j+d_count,k+j+d_count + 1
  !          T(k+j+d_count) = 0.5*X(j)
  !          d_count = d_count + 1
  !          T(k+j+d_count) = 0.5*(X(j) + X(j+1))
  !          d_count = d_count + 1
  !          T(k+j+d_count) = 0.5*(X(j) + 1)
  !       else
  !          print *,'should not be here'
  !    !      T(k+j+d_count) = X(j+1)
  !       end if

  !    end do
  !    print *, 'caled knots are:',T
  ! else
  !    T(k:Nctl) = temp(k:Nctl) ! Copy the temp vector verbatem
  ! end if
  T(4) = X(2)/2
  T(5) = X(2)
  T(6) = 0.5*(X(2)+1)
  !T(5) = X(2)/2
  !T(6) = (X(2)+1)/2
  print *,'knots are :',T
  print *,'X is :',X
end subroutine knots_interp
