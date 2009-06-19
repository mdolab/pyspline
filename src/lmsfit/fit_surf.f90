subroutine fit_surf(Nsurf,Nu,Nv,Nctlu,Nctlv,Ncon,A,X,B,D,ctl)
!-----------------------------------------------------------------------
! Purpose: Fit_surf performs a constrained Least Mean Squares Fit on a
! set of x-y-z points in 3-space defined by a parametric values u,v. 
!
! Inputs:
!
! A: This is the matrix of partial derivaties of size M by N
!
! Nu: This is the number of datapoints in U we are fitting
! Nv: This is the number of datapoints in V we are fitting
!
! N: This is the number of control points we have
! 

! Parameters

! Input/Output Variables
implicit none
integer, intent(in) :: Nsurf ! Number of surfaces
integer, intent(in) :: Nu    ! Number of points in u
integer, intent(in) :: Nv    ! Number of points in v
integer, intent(in) :: Nctlu ! Number of control points in u
integer, intent(in) :: Nctlv ! Number of control points in v
integer, intent(in) :: Ncon  ! Number of LINEAR constraints

double precision, intent(in) :: A(Nsurf,Nu*Nv,Nctlu*Nctlv)
double precision, intent(in) :: X(Nsurf,Nu,Nv,3)
double precision, intent(in) :: B(Ncon,Nsurf*3*Nctlu*Nctlv)
double precision, intent(in) :: D(Ncon)

double precision, intent(out):: ctl(Nsurf,Nctlu,Nctlv,3)

! Local Variables
double precision :: x_temp(Nu*Nv,1)
integer isurf,idim,i,j

!-- variables needed for QR decomposition
integer :: lwork, info
double precision :: work(2*Nctlu*Nctlv*Nu*Nv)
double precision err
character trans

external dgglse ! Constrained
external dgels ! Unconstrained

! Code for doing unconstrained least squares
!print *,'In fit_surf'
do isurf=1,Nsurf
   do idim = 1,3
      do i=1,Nu
         do j=1,Nv
            X_temp((i-1)*Nv + j,1) = X(isurf,i,j,idim)
         end do
      end do
      !print *,'isurf,idim',isurf,idim
      !print *,X_temp
      trans = 'n'
      lwork = 2*Nctlu*Nctlv*Nu*Nv
      !call dgels( TRANS, M,     N,      NRHS, A,          LDA,  B,    LDB, WORK, LWORK, INFO )
      call dgels(trans  ,Nu*Nv,Nctlu*Nctlv,1,  A(isurf,:,:),Nu*Nv, X_temp,Nu*Nv,    work, lwork,  info)
      do i =1,Nctlu
         do j = 1,Nctlv
            ctl(isurf,i,j,idim) = X_temp((i-1)*Nctlv + j,1)
         end do
      end do
   end do
end do

!print *,'done fit_surf'
end subroutine fit_surf
