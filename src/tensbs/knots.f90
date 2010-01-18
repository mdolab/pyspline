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

  d = dble(N/(Nctl-k+1.0))

  do j=1,Nctl-k+1
     i = floor(j*d)
     alpha = j*d-i
     T(k+j) = (1-alpha)*X(i) + alpha*X(i+1)
  end do

end subroutine knots



!   if (Nctl == N) then
!      call bknot(X,N,K,T)
!      return
!   end if

!   ! We now have Nctl-K knots left to distribute
!   ! We wich to distribute them according to positions of data
!   ! 
!   Nremain = Nctl-K

!   if (Nremain == 0) then
!      return
!   else
     
!      if (mod(Nremain,2) == 0) then
!         !print *,'even remaining'

!         istart = (Nctl+k)/2 

!         frac = 1.0/(Nremain + 1)
!         do i =1,Nremain/2
!            ! Lower
!            temp = frac*i*N
!            index = int(floor(temp))
!            s = temp-index

!            T(k+i) = X(index)+ s*(X(index+1)-X(index))

!            ! Upper
!            temp =(1-frac*i)*N+1
!            index = int(floor(temp))
!            s = temp-index

!            T(Nctl-i+1) = X(index) + s*(X(index+1)-X(index))
!         end do
!      else
!         !print *,'odd remaining'
        
!         istart= (Nctl+k+1)/2
!         frac = 1.0/(Nremain + 1)

!         ! Do The middle...always at FRAC=0.5
!         temp = 0.5*N+0.5
!         index = int(floor(temp))
!         s = temp-index
!         T(istart) = X(temp) + s*(X(temp+1)-X(temp))

!         do i =1,(Nremain-1)/2

!            ! Lower
!            temp = (0.5-frac*i)*N
!            index = int(floor(temp))
!            s = temp-index

!            T(istart-i) = X(index)+ s*(X(index+1)-X(index))

!            ! Upper
!            temp = (0.5+frac*i)*N
!            index = int(floor(temp))
!            s = temp-index

!            T(istart+i) = X(index)+ s*(X(index+1)-X(index))

!         end do
!      end if
!   end if
! end subroutine knots
