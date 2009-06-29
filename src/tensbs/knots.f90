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
  double precision frac,s
  INTEGER  I, J, Nremain, istart, highI,lowI,mid,temp
  ! ----------------------------
  !  PUT K KNOTS AT EACH ENDPOINT -- Knot a knot conditions
  !  ----------------------------
  !
  print *,'in knots'
  DO  J=1,K
     T(J) = X(1) ! Left
     T(Nctl+J) = X(N) ! right
  end do

  ! We now have Nctl-K knots left to distribute
  ! We wich to distribute them according to positions of data
  ! 
  Nremain = Nctl-K
print *,'N:',N
  if (Nremain == 0) then
     return
  else
     
     if (mod(Nremain,2) == 0) then
        print *,'even remaining'

        istart = (Nctl+k)/2

        frac = 1.0/(Nremain + 1)

        do i =1,Nremain/2
           
           temp = int(floor(frac*i*N))
           s = frac*i*N-temp
           print *,'index:',k+i,temp,s
!  !          print *,'temp,s',temp,s
           T(k+i) = X(temp) + s*(X(temp+1)-X(temp))
!            print *,'X(temp);',X(temp)
           ! figure out where to put values
           
        end do
        do i=1,Nremain/2
           temp =int(floor((1-frac*i)*N))+1
           s =  (1-frac*i)*N+1-temp
           print *,'index:',Nctl-i+1,temp,s
! !           print *,'temp,s:',temp,s
           T(Nctl-i+1) = X(temp) + s*(X(temp+1)-X(temp))
!            print *,'X(temp)',X(temp)


        end do
     else
        print *,'odd remaining'
     end if
  end if












! ! Even number of data points
!            istart = (Nctl+k)/2 !should be whole
!            mid = N/2
!            frac = 1.0/(Nctl-k+2)
!            print *,'frac:',frac
!            do i=1,Nremain/2
!               ! Must set two vals here
!               ! Frac tells us what how far to go to get values
!               highI = mid + int(N/2)*frac
!               lowI  = mid - int(N/2)*frac
              
!               print *,'lowI,highI',lowI,highI,i,istart
!               print *,'istart-1+1',istart-1+i
!               print *,'istart+i',istart+i
!               T(istart-i+1) = X(lowI)
!               T(istart+i  ) = X(highI)

!            end do
           
    

end subroutine knots
