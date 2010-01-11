subroutine compute_curve(s,X,t,k,n,nctl,ndim,coef)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: compute_curve is the main function for the generation fo NURBS curves
  !               It does both interpolating and LMS fits, variable NURBS weights as well as
  !               Hoschek's Parameter Correction.
  !
  !     Description of Arguments
  !     Input
  !     s       - Real, Vector of s coordinates, length n
  !     X       - Real, Array of X values to fit, Size (n,ndim)
  !     t       - Real,Knot vector. Length nctl+k
  !     k       - Integer,order of spline
  !     nctl    - Integer,Number of control points
  !     ndim    - Integer, spatial dimension of curve
  !
  !     Ouput 
  !     coef    - Real,Array of NURBS coefficients and weights. Size (nctl,ndim+1)  

  implicit none
  
  ! Input
  integer         , intent(in)          :: k,nctl,ndim,n
  double precision, intent(in)          :: X(n,ndim)
  double precision, intent(inout)       :: s(n)
  double precision, intent(in)          :: t(nctl+k)

  ! Output
  double precision, intent(out)         :: coef(nctl,ndim+1)

  ! Working
  integer                               :: i,idim,inbv
  integer                               :: niter,iter
  double precision                      :: work(3*k)
  double precision                      :: weight
!  double precision                      :: Jac(n,nctl*(ndim+1))
  double precision                      :: Xcopy(n,ndim)
  double precision                      :: Jac(n,nctl)
  integer                               :: info,lwork
  double precision                      :: work2(nctl+max(nctl,ndim))
  double precision                      :: length 

  double precision                      :: poly_length
  
  lwork = nctl + max(nctl,ndim)

  print *,'Welcome to Compute Curve'

  niter = 1
  coef(:,ndim+1) = 1.0 ! Initialize all the weights to 1
  Jac(:,:) = 0.0

  ! First we want to do a Linear Least Squares to get a good starting point

  Xcopy(:,:) = X(:,:)   ! Copy 'X' to Xcopy
  call curve_jacobian(Jac,n,nctl,ndim,t,k,s,coef)
  call DGELS('N', n, nctl, ndim, Jac, n, Xcopy, n, WORK, LWORK, INFO )
  do i=1,nctl
     do idim=1,ndim
        coef(i,idim) = Xcopy(i,idim)
     end do
  end do




  length = poly_length(X,n,ndim)
  print *,'length:',length

  coef(:,ndim+1) = 1.0 ! Set all the weights to 1
  
end subroutine compute_curve


subroutine curve_jacobian(Jac,n,nctl,ndim,t,k,s,coef)

  implicit none
  integer                               :: k,nctl,ndim,n
  double precision                      :: Jac(n,nctl)
  double precision                      :: s(n)
  double precision                      :: t(nctl+k)
  double precision                      :: coef(nctl,ndim+1)
  double precision                      :: vnikx(k),work(2*k)
  integer                               :: i,j,idim
  integer                               :: ilo,ileft,mflag,iwork
  
  ilo = 1
  do i=1,n
     call intrv(t,nctl+k,s(i),ilo,ileft,mflag)
     if (mflag == 0) then
        call bspvn(t,k,k,1,s(i),ileft,vnikx,work,iwork)
     else if (mflag == 1) then
        ileft = nctl
        vnikx(:) = 0.0
        vnikx(k) = 1.0
     end if

     do j=1,k
        Jac(i,ileft-k+j) = vnikx(j)
     end do
  end do

end subroutine curve_jacobian

function poly_length(X,n,ndim)
  ! Compute the length of the spatial polygon
  implicit none

  double precision                :: X(n,ndim)
  integer                         :: n,ndim,i,idim
  double precision                :: dist
  double precision poly_length
  poly_length = 0.0
  do i=1,n-1
     dist = 0.0
     do idim=1,ndim
        dist = dist + (X(i,idim)-X(i+1,idim))**2
     end do
     poly_length = poly_length + sqrt(dist)
  end do
  
end function poly_length
