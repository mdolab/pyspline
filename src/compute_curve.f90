! subroutine compute_curve(s,X,t,k,n,nctl,ndim,coef,niter,tol,rms)

!   !***DESCRIPTION
!   !
!   !     Written by Gaetan Kenway
!   !
!   !     Abstract: compute_curve is the main function for the
!   !               generation fo B-spline curves It does both
!   !               interpolating and LMS fits, as well as Hoschek's 
!   !               Parameter Correction.
!   !
!   !     Description of Arguments
!   !     Input
!   !     s       - Real, Vector of s coordinates, length n
!   !     X       - Real, Array of X values to fit, Size (n,ndim)
!   !     t       - Real,Knot vector. Length nctl+k
!   !     k       - Integer,order of spline
!   !     nctl    - Integer,Number of control points
!   !     ndim    - Integer, spatial dimension of curve
!   !     niter   - Integer, number of hoscek para-correction iterations
!   !     tol     - Integer, relative tolerance for convergence
!   !
!   !     Ouput coef - Real,Array of B-spline coefficients. Size (nctl,ndim)

!   use lms_jacobian
!   use  lsqrModule,        only : LSQR
!   use  lsqrCheckModule,   only : Acheck, xcheck

!   use precision
!   implicit none
  
!   ! Input
!   integer         , intent(in)          :: k,nctl,ndim,n
!   real(kind=realType), intent(in)          :: X(n,ndim)
!   real(kind=realType), intent(inout)       :: s(n)
!   real(kind=realType), intent(in)          :: t(nctl+k)
!   integer         , intent(in)          :: niter
!   real(kind=realType), intent(in)          :: tol

!   ! Output
!   real(kind=realType), intent(inout)       :: coef(nctl,ndim)
!   real(kind=realType), intent(out)         :: rms
!   ! Working
!   integer                               :: i,idim,iter
!   real(kind=realType)                      :: weight
!   real(kind=realType)                      :: length,res,norm,tot,val(ndim)
!   integer                               :: istop,itn
!   real(kind=realType)                      :: Anorm,Acond,rnorm, Arnorm,xnorm
!   ! Functions called
!   real(kind=realType)                      :: poly_length,floor,compute_rms_curve

!   !Initialization
!   length = poly_length(X,n,ndim)  
!   call setup_jacobian(n,nctl,k)
!   call curve_jacobian_linear(t,k,s,n,nctl)

!   do iter=1,niter
!      ! Solve
!      idim = 1
!      do idim=1,ndim
!         call LSQR( n, Nctl, Aprod1, Aprod2,X(:,idim),0.0, .False., &
!              coef(:,idim), vals, 1e-12, 1e-12, 1e8, Nctl*4,0,&
!              istop, itn, Anorm, Acond, rnorm, Arnorm, xnorm )
!      end do
!      call curve_para_corr(t,k,s,coef,nctl,ndim,length,n,X)
!      rms = compute_rms_curve(t,k,s,coef,nctl,ndim,length,n,X)
!      call curve_jacobian_linear(t,k,s,n,nctl)
!      ! Do convergence Check to break early
!   end do
 
!   call kill_jacobian()

! end subroutine compute_curve

! subroutine curve_jacobian_linear(t,k,s,n,nctl)
!   use lms_jacobian

!   use precision
!   implicit none
!   ! Input
!   real(kind=realType)     ,intent(in)      :: t(nctl+k)
!   real(kind=realType)     ,intent(in)      :: s(n)
!   integer              ,intent(in)      :: k,nctl,n
!   ! Output -> in lms_jacobian module

!   real(kind=realType)                      :: vnikx(k),work(2*k)
!   integer                               :: i,j,counter
!   integer                               :: ilo,ileft,mflag,iwork
  
!   ilo = 1
!   counter = 1
!   do i=1,n
!      call intrv(t,nctl+k,s(i),ilo,ileft,mflag)
!      if (mflag == 0) then
!         call bspvn(t,k,k,1,s(i),ileft,vnikx,work,iwork)
!      else if (mflag == 1) then
!         ileft = nctl
!         vnikx(:) = 0.0
!         vnikx(k) = 1.0
!      end if

!      row_ptr(i) = counter
!      do j=1,k
!         col_ind(counter) = ileft-k+j
!         vals(counter) =vnikx(j)
!         counter = counter + 1
!      end do
!   end do
!   row_ptr(n+1) = counter
! end subroutine curve_jacobian_linear

subroutine curve_jacobian_wrap(s,sd,t,k,nctl,n,nd,vals,row_ptr,col_ind)
  use precision
  implicit none
  ! Input
  real(kind=realType)     ,intent(in)      :: t(nctl+k),s(n),sd(nd)
  integer              ,intent(in)      :: k,nctl,n,nd
  ! Output
  real(kind=realType)     ,intent(inout)   :: vals((n+nd)*k)
  integer              ,intent(inout)   :: row_ptr(n+nd+1)
  integer              ,intent(inout)   :: col_ind((n+nd)*k)
  real(kind=realType)                      :: basisu(k),work((k+1)*(k+2)/2),vdikx(k,2)
  integer                               :: i,j,counter
  integer                               :: ilo,ileft,mflag,iwork
  ilo = 1
  counter = 1
  do i=1,n ! Do the values first 
     call intrv(t,nctl+k,s(i),ilo,ileft,mflag)
     if (mflag == 1) then
        ileft = ileft - k
     end if
     call basis(t,nctl,k,s(i),ileft,basisu)

     row_ptr(i) = counter-1
     do j=1,k
        col_ind(counter) = ileft-k+j-1
        vals(counter) = basisu(j)
        counter = counter + 1
     end do
  end do
  do i=1,nd ! Do the derivatives next
     call intrv(t,nctl+k,sd(i),ilo,ileft,mflag)
     if (mflag == 1) then
        ileft = ileft - k
     end if
     call bspvd(t,k,2,sd(i),ileft,k,vdikx,work)
     row_ptr(i+n) = counter-1
     do j=1,k
        col_ind(counter) = ileft-k+j-1
        vals(counter) =vdikx(j,2)
        counter = counter + 1
     end do
  end do
  row_ptr(n+nd+1) = counter-1
end subroutine curve_jacobian_wrap

subroutine constr_jac(A_val,A_row_ptr,A_col_ind,B_val,B_row_ptr,B_col_ind,C_val,C_row_ptr,C_col_ind, &
     Am,An,Cm,Annz,Bnnz,Cnnz,J_val,J_col_ind,J_row_ptr)
  use precision
  implicit none
  ! Input
  integer         , intent(in)   :: Am,An,Cm,Annz,Bnnz,Cnnz
  real(kind=realType), intent(in)   :: A_val(Annz),B_val(Bnnz),C_val(Cnnz)
  integer         , intent(in)   :: A_col_ind(Annz),B_col_ind(Bnnz),C_col_ind(Cnnz)
  integer         , intent(in)   :: A_row_ptr(Am+1),B_row_ptr(Am+1),C_row_ptr(Cm+1)

  ! Output
  real(kind=realType), intent(out)  :: J_val(Annz+Bnnz+Cnnz)
  integer         , intent(out)  :: J_col_ind(Annz+Bnnz+Cnnz)
  integer         , intent(out)  :: J_row_ptr(Am+Cm+1)

  ! Local 
  integer                        :: i,j,counter
  ! This functions assembes the following CSR matrix:
  ! J = [A    B]
  !     [C    0]

  ! Now assmeble the full jacobain
  counter = 1
  J_row_ptr(1) = 0
  do i =1,Am
     do j=1,A_row_ptr(i+1)-A_row_ptr(i)
        J_val(counter) =     A_val(A_row_ptr(i)+j)
        J_col_ind(counter) = A_col_ind(A_row_ptr(i)+j)
        counter = counter + 1
     end do
     do j=1,B_row_ptr(i+1)-B_row_ptr(i)
        J_val(counter) = B_val(B_row_ptr(i)+j)
        J_col_ind(counter) = B_col_ind(B_row_ptr(i)+j) + An
        counter = counter + 1
     end do
     J_row_ptr(i+1) = counter - 1
  end do
  do i =1,Cm
     do j=1,C_row_ptr(i+1)-C_row_ptr(i)
        J_val(counter) = C_val(C_row_ptr(i)+j)
        J_col_ind(counter) = C_col_ind(C_row_ptr(i)+j)
        counter = counter + 1
     end do
     J_row_ptr(i+1+am) = counter - 1
  end do
end subroutine constr_jac

function poly_length(X,n,ndim)
  ! Compute the length of the spatial polygon
  use precision
  implicit none
  !Input
  real(kind=realType) ,intent(in)    :: X(ndim,n)
  integer          ,intent(in)    :: n,ndim
  integer                         :: i,idim
  real(kind=realType)                :: dist
  real(kind=realType) poly_length
  poly_length = 0.0
  do i=1,n-1
     dist = 0.0
     do idim=1,ndim
        dist = dist + (X(idim,i)-X(idim,i+1))**2
     end do
     poly_length = poly_length + sqrt(dist)
  end do
  
end function poly_length

subroutine curve_para_corr(t,k,s,coef,nctl,ndim,length,n,X)

  ! Do Hoschek parameter correction
  use precision
  implicit none
  ! Input/Output
  integer           ,intent(in)      :: k,nctl,ndim,n
  real(kind=realType)  ,intent(in)      :: t(nctl+k)
  real(kind=realType)  ,intent(inout)   :: s(n)
  real(kind=realType)  ,intent(in)      :: coef(ndim,nctl)
  real(kind=realType)  ,intent(in)      :: X(ndim,n)
  real(kind=realType)  ,intent(in)      :: length
  ! Working
  integer                            :: i,j,max_inner_iter
  real(kind=realType)                   :: D(ndim),D2(ndim)
  real(kind=realType)                   :: val(ndim),deriv(ndim)
  real(kind=realType)                   :: c,s_tilde,norm

  max_inner_iter = 10
  do i=2,n-1
     call eval_curve(s(i),t,k,coef,nctl,ndim,1,val)
     call eval_curve_deriv(s(i),t,k,coef,nctl,ndim,deriv)
          
     D = X(:,i)-val
     c = dot_product(D,deriv)/norm(deriv,ndim)

     inner_loop: do j=1,max_inner_iter

        s_tilde = s(i)+ c*(t(nctl+k)-t(1))/length
        call eval_curve(s_tilde,t,k,coef,nctl,ndim,1,val)
        D2 = X(:,i)-val
        if (norm(D,ndim) .ge. norm(D2,ndim)) then
           s(i) = s_tilde
           exit inner_loop
        else
           c = c*0.5
        end if
     end do inner_loop
  end do

end subroutine curve_para_corr

function norm(X,n)
  ! Compute the L2 nomr of X
  use precision
  implicit none
  integer, intent(in) :: n
  real(kind=realType),intent(in) :: X(n)

  real(kind=realType) :: norm

  integer            :: i
  norm = 0.0
  do i=1,n
     norm = norm + X(i)**2
  end do
  norm = sqrt(norm)
end function norm

function compute_rms_curve(t,k,s,coef,nctl,ndim,length,n,X)
  ! Compute the rms
  use precision
  implicit none
  ! Input/Output
  real(kind=realType)  ,intent(in)      :: t(k+nctl)
  real(kind=realType)  ,intent(in)      :: s(n)
  real(kind=realType)  ,intent(in)      :: coef(ndim,nctl)
  integer           ,intent(in)      :: k,nctl,ndim,n
  real(kind=realType)  ,intent(in)      :: X(ndim,n)
  real(kind=realType)  ,intent(in)      :: length
  real(kind=realType)                   :: compute_rms_curve 

  ! Working
  integer                            :: i,idim
  real(kind=realType)                :: val(ndim),D(ndim)
  
  compute_rms_curve = 0.0
  do i=1,n
     call eval_curve(s(i),t,k,coef,nctl,ndim,1,val)
     D = val-X(:,i)
     do idim=1,ndim
        compute_rms_curve = compute_rms_curve + D(idim)**2
     end do
  end do
  compute_rms_curve = sqrt(compute_rms_curve/n)
end function compute_rms_curve
