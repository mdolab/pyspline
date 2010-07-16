subroutine insertKnot(u,r,t,k,coef,nctl,ndim,t_new,coef_new,ileft)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract insertKnot inserts a knot u into the curve, s times
  !
  !     Description of Arguments
  !     Input
  !     u       - knot to insert
  !     r       - Insert r times
  !     t       - Real,Knot vector. Length nctl+k
  !     k       - Integer,order of B-spline
  !     coef    - Real,Array of B-spline coefficients  Size (ndim,nctl)
  !     nctl   - Integer,Number of control points
  !
  !     Ouput 
  !     t_new - Real, vector of lenght(nctl+k+r)
  !     coef_new - Array of new cofficients size(ndim,nctl+r)
  implicit none
  ! Input
  integer         , intent(inout)       :: r
  integer         , intent(in)          :: k,nctl,ndim
  double precision, intent(in)          :: u
  double precision, intent(in)          :: t(nctl+k)
  double precision, intent(in)          :: coef(ndim,nctl)
  ! Output
  double precision, intent(out)         :: t_new(nctl+k+r)
  double precision, intent(out)         :: coef_new(ndim,nctl+r)
  integer         , intent(out)         :: ileft
  ! Working
  integer                               :: mflag,ilo,s,i,j,L
  double precision                      :: alpha,temp(ndim,k)

  ilo = 1
  call intrv(t,nctl+k,u,ilo,ileft,mflag)
  ! Compute its multiplicity
  s = 0 ! Knot multiplicity
  mult: do i=0,k-1
     if (t(ileft-i) == u) then 
        s = s + 1
     else
        exit mult
     end if
  end do mult
  
  ! We need to make sure that the requested multipliity r, plus 
  ! the actual multiplicity of this know is less than (k-1)

  if (s + r + 1 > k) then
     r = k-1-s
  end if
  ! --------- Load New Knot Vector -------
  do i=1,ileft
     t_new(i) = t(i)
  end do

  do i=1,r
     t_new(ileft+i) = u
  end do

  do i=ileft+1,nctl+k
     t_new(i+r) = t(i)
  end do

  ! -------- Save unaltered Control Points

  do i=1,ileft-(k-1)
     coef_new(:,i) = coef(:,i)
  end do
  
  do i=ileft,nctl
     coef_new(:,i+r) = coef(:,i)
  end do
  
  do i=1,k-s
     temp(:,i) = coef(:,ileft-k+i)
  end do

  do j=1,r
     L = ileft-k+j+1
     do i=0,k-1-j-s
        alpha = (u-t(L+i))/(t(i+ileft+1)-t(L+i))
        temp(i+1,:) = alpha*temp(:,i+2) + (1.0-alpha)*temp(:,i+1)
     end do
     coef_new(:,L) = temp(:,1)
     coef_new(:,ileft+r-j) = temp(:,k-j)
  end do
  do i=L+1,ileft-1
     coef_new(:,i) = temp(:,i-L+1)
  end do
end subroutine insertKnot
 
