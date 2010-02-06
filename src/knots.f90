subroutine knots(X,N,Nctl,K,T)
  !  --------------------------------------------------------------------
  !  knots chooses a knot sequence,T, for data,X, of length N, with Nctl
  !  control points, and a spline of order k 
  !  --------------------------------------------------------------------

  implicit none
  integer, intent(in)     :: N
  integer, intent(in)     :: Nctl
  integer, intent(in)     :: K

  double precision, intent(in)  :: X(N)
  double precision, intent(out) :: T(Nctl+K)
  !
  ! LOCAL VARIABLES
  !
  double precision d,alpha
  INTEGER  I, J
  ! ----------------------------
  !  PUT K KNOTS AT EACH ENDPOINT -- Knot a knot conditions
  !  ----------------------------
  !

  DO  J=1,K
     T(J) = X(1) ! Left
     T(Nctl+J) = X(N) ! right
  end do

  if (N .ne. Nctl) then

     d = dble(N/(Nctl-k+1.0))

     do j=1,Nctl-k+1
        i = floor(j*d)
        alpha = j*d-i
        T(k+j) = (1-alpha)*X(i) + alpha*X(i+1)
     end do
   else  ! N == Nctl
      do j = 1,N-k+1
         T(k+j) = X(j+1)
      end do
   end if
      


end subroutine knots
