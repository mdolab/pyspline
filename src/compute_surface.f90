! subroutine compute_surface(X,u,v,tu,tv,ku,kv,nctlu,nctlv,nu,nv,ndim,coef,niter,tol,rms)
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
!   !     X       - Real, Array of X values to fit, Size (nu,nv,ndim)  
!   !     u       - Real, Array of u coordinates, Size nu x nv
!   !     v       - Real, Array of v coordinates, Size nu x nv
!   !     tu      - Real,Knot vector in u. Length nctlu+ku
!   !     tv      - Real,Knot vector in v. Length nctlv+kv
!   !     ku      - Integer, order of B-spline in u
!   !     kv      - Integer, order of B-spline in v
!   !     nctlu   - Integer,Number of control points in u
!   !     nctlv   - Integer,Number of control points in v
!   !     nu      - Integer, Number of data points in u
!   !     nv      - Integer, Number of data points in v
!   !     ndim    - Integer, spatial dimensions
!   !     coef    - Array, B-spline coefficients. Size nctlu x nctlv x ndim
!   !     niter   - Integer, number of hoscek para-correction iterations
!   !     tol     - Integer, relative tolerance for convergence
!   !
!   !     Ouput 
!   !     rms     - Real, final RMS value

!   use lms_jacobian
!   use lsqrModule,        only : LSQR
!   use lsqrCheckModule,   only : Acheck, xcheck

!   use precision
!   implicit none
  
!   ! Input/Output
!   real(kind=realType), intent(in)          :: X(nu,nv,ndim)
!   integer         , intent(in)          :: ku,kv,nctlu,nctlv,nu,nv,ndim
!   real(kind=realType), intent(inout)       :: u(nu,nv),v(nu,nv)
!   real(kind=realType), intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
!   real(kind=realType), intent(inout)       :: coef(nctlu,nctlv,ndim)
!   integer         , intent(in)          :: niter
!   real(kind=realType), intent(in)          :: tol
!   real(kind=realType), intent(out)         :: rms

!   ! Working
!   integer                               :: i,idim,iter
!   integer                               :: istop,itn
!   real(kind=realType)                      :: Anorm,Acond,rnorm, Arnorm,xnorm
!   ! Functions called
!   real(kind=realType)                      :: poly_length,floor,compute_rms_surface

!   !Initialization
!   call setup_jacobian(nu*nv,nctlu*nctlv,ku*kv)
!   call surface_jacobian_linear(u,v,tu,tv,ku,kv,nctlu,nctlv,nu,nv)
!   rms = 0.0


!   do iter=1,niter
! !      ! Solve
!       idim = 1
!       do idim=1,ndim
!          call LSQR( nu*nv, Nctlu*Nctlv, Aprod1, Aprod2,X(:,:,idim),.0, .False., &
!               coef(:,:,idim), vals, 1e-12, 1e-12, 1e8, Nctlu*Nctlv*10,0,&
!               istop, itn, Anorm, Acond, rnorm, Arnorm, xnorm )
!          print *,'istop,itn:',istop,itn
!          print *,'Anorm,Acond:',Anorm,Acond,Arnorm

!       end do
!       !call surface_para_corr(t,k,s,coef,nctl,ndim,length,n,X)
!       rms = compute_rms_surface(tu,tv,ku,kv,u,v,coef,nctlu,nctlv,ndim,nu,nv,X)
!       call surface_jacobian_linear(u,v,tu,tv,ku,kv,nctlu,nctlv,nu,nv)
!      ! Do convergence Check to break early

!   end do
!   call kill_jacobian()

! end subroutine compute_surface

! subroutine surface_jacobian_linear(u,v,tu,tv,ku,kv,nctlu,nctlv,nu,nv)

!   !***DESCRIPTION
!   !
!   !     Written by Gaetan Kenway
!   !
!   !     Abstract surface_jacobian_linear computes jacobian in rcs storage
!   !
!   !     Description of Arguments
!   !     Input
!   !     u       - Real, u coordinates, size: nu x nv
!   !     v       - Real, v coordinates, size: nu x nv
!   !     tu      - Real,Knot vector in u. Length nctlu+ku
!   !     tv      - Real,Knot vector in v. Length nctlv+kv
!   !     ku      - Integer, order of B-spline in u
!   !     kv      - Integer, order of B-spline in v
!   !     nctlu   - Integer,Number of control points in u
!   !     nctlv   - Integer,Number of control points in v
!   !     nu      - Integer, Number of data points in u
!   !     nv      - Integer, Number of data points in v
!   !
!   !     Ouput -> lms_jacobian module

!   use lms_jacobian

!   use precision
!   implicit none
!   ! Input
!   integer         , intent(in)          :: ku,kv,nctlu,nctlv,nu,nv
!   real(kind=realType), intent(in)          :: u(nu,nv),v(nu,nv)
!   real(kind=realType), intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)

!   ! Working
!   real(kind=realType)                      :: vniku(ku),worku(2*ku)
!   integer                               :: ilou,ileftu,mflagu

!   real(kind=realType)                      :: vnikv(kv),workv(2*kv)
!   integer                               :: ilov,ileftv,mflagv

!   integer                               :: i,j,ii,jj,counter,iwork

!   ilou = 1
!   ilov = 1
!   counter = 1
!   do i=1,nu
!      do j = 1,nv
!         ! Get u interval
!         call intrv(tu,nctlu+ku,u(i,j),ilou,ileftu,mflagu)
!         if (mflagu == 0) then
!            call bspvn(tu,ku,ku,1,u(i,j),ileftu,vniku,worku,iwork)
!         else if (mflagu == 1) then
!            ileftu = nctlu
!            vniku(:) = 0.0
!            vniku(ku) = 1.0
!         end if

!         ! Get v interval
!         call intrv(tv,nctlv+kv,v(i,j),ilov,ileftv,mflagv)
!         if (mflagv == 0) then
!            call bspvn(tv,kv,kv,1,v(i,j),ileftv,vnikv,workv,iwork)
!         else if (mflagv == 1) then
!            ileftv = nctlv
!            vnikv(:) = 0.0
!            vnikv(kv) = 1.0
!         end if
        
!         row_ptr( (i-1)*nv + j ) = counter
!         do ii=1,ku
!            do jj = 1,kv
!               !rows(counter) = (i-1)*nv + j -1
!               ! cols(counter) = (ileftu-ku+ii-1)*Nctlv + (ileftv-kv+jj-1)
!               col_ind(counter) = (ileftu-ku+ii-1)*Nctlv + (ileftv-kv+jj)
!               vals(counter) = vniku(ii)*vnikv(jj)
!               counter = counter + 1
!            end do
!         end do
!      end do
!   end do
!   row_ptr(nu*nv+1) = counter
! end subroutine surface_jacobian_linear


subroutine surface_jacobian_wrap(u,v,tu,tv,ku,kv,nctlu,nctlv,nu,nv,vals,row_ptr,col_ind)

  use precision
  implicit none
  ! Input
  integer         , intent(in)          :: ku,kv,nctlu,nctlv,nu,nv
  real(kind=realType), intent(in)          :: u(nv,nu),v(nv,nu)
  real(kind=realType), intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
  real(kind=realType), intent(out)         :: vals(nu*nv*ku*kv)
  integer         , intent(out)         :: col_ind(nu*nv*ku*kv),row_ptr(nu*nv+1)
  ! Working
  real(kind=realType)                      :: basisu(ku),basisv(kv)
  integer                               :: ilou,ileftu,mflagu
  integer                               :: ilov,ileftv,mflagv

  integer                               :: i,j,ii,jj,iwork,counter

  ilou = 1
  ilov = 1
  counter = 1
  !print *,'uknot:',tu
  !print *,'vknot:',tv


  do i=1,nu
     do j = 1,nv
        ! Get u interval
        call intrv(tu,nctlu+ku,u(j,i),ilou,ileftu,mflagu)
        if (mflagu == 0) then
           call basis(tu,nctlu,ku,u(j,i),ileftu,basisu)
        else if (mflagu == 1) then
           ileftu = nctlu
           basisu(:) = 0.0
           basisu(ku) = 1.0
        end if
        
        ! Get v interval
        call intrv(tv,nctlv+kv,v(j,i),ilov,ileftv,mflagv)
        if (mflagv == 0) then
           call basis(tv,nctlv,kv,v(j,i),ileftv,basisv)
        else if (mflagv == 1) then
           ileftv = nctlv
           basisv(:) = 0.0
           basisv(kv) = 1.0
        end if

        row_ptr( (i-1)*nv + j  ) = counter-1
        do ii=1,ku
           do jj = 1,kv
              col_ind(counter) = (ileftu-ku+ii-1)*Nctlv + (ileftv-kv+jj-1)
              vals(counter) = basisu(ii)*basisv(jj)
              counter = counter + 1
              
           end do
        end do
     end do
  end do
  row_ptr(nu*nv+1) = counter-1
  
end subroutine surface_jacobian_wrap

subroutine surface_para_corr(tu,tv,ku,kv,u,v,coef,nctlu,nctlv,ndim,nu,nv,X,rms)

  ! Do Hoschek parameter correction
  use precision
  implicit none
  ! Input/Output
  real(kind=realType)  ,intent(in)      :: tu(ku+nctlu),tv(kv+nctlv)
  real(kind=realType)  ,intent(inout)   :: u(nv,nu),v(nv,nu)
  real(kind=realType)  ,intent(in)      :: coef(ndim,nctlv,nctlu)
  integer           ,intent(in)      :: ku,kv,nctlu,nctlv,ndim,nu,nv
  real(kind=realType)  ,intent(in)      :: X(ndim,nv,nu)
  real(kind=realType)  ,intent(out)     :: rms
  ! Working
  integer                            :: i,j,jj,max_inner_iter
  real(kind=realType)                   :: lengthu,lengthv
  real(kind=realType)                   :: D(ndim),D2(ndim)
  real(kind=realType)                   :: val(ndim),deriv(ndim,2),deriv2(ndim,2,2)
  real(kind=realType)                   :: delta_c,delta_d,u_tilde,v_tilde
  integer                            :: adj_u,adj_v
  real(kind=realType)                   :: A(2,2),ki(2),delta(2)
  !Functions
  real(kind=realType)                   :: norm,poly_length

  max_inner_iter = 10
  rms = 0.0


  do i=2,nu-1
     do j = 2,nv-2
        call eval_surface(u(j,i),v(j,i),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
        call eval_surface_deriv(u(j,i),v(j,i),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv)
        call eval_surface_deriv2(u(j,i),v(j,i),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv2)

        D = val-X(:,j,i)

        A(1,1) = norm(deriv(:,i),ndim)**2 + dot_product(D,deriv2(:,1,1))
        A(1,2) = dot_product(deriv(:,1),deriv(:,2)) + dot_product(D,deriv2(:,1,2))
        A(2,1) = A(1,2)
        A(2,2) = norm(deriv(:,2),ndim)**2 + dot_product(D,deriv2(:,2,2))
        
        ki(1) = -dot_product(D,deriv(:,1))
        ki(2) = -dot_product(D,deriv(:,2))
        
        call solve_2by2(A,ki,delta)
        
        if (j .eq.1 .or. j .eq. nv) then
           delta(1) = 0.0
        end if
        if (i .eq.1 .or. i .eq. nu) then
           delta(2) = 0.0
        end if
        inner_loop: do jj=1,max_inner_iter
           u_tilde = u(j,i) + delta(1)
           v_tilde = v(j,i) + delta(2)

           call eval_surface(u_tilde,v_tilde,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
           D2 = val-X(i,j,:)
           if (norm(D,ndim) .ge. norm(D2,ndim)) then
              u(j,i) = u_tilde
              v(j,i) = v_tilde
              exit inner_loop
           else
              delta = delta*0.5
           end if
        end do inner_loop
     end do
  end do

  ! Lets redo the full RMS
  rms = 0.0


  do i=1,nu
     do j=1,nv
        call eval_surface(u(j,i),v(j,i),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
        D = X(i,j,:)-val
        rms = rms + dot_product(D,D)
     end do
  end do
  rms = sqrt(rms/(nu*nv))
  

end subroutine surface_para_corr

function compute_rms_surface(tu,tv,ku,kv,u,v,coef,nctlu,nctlv,ndim,nu,nv,X)
 ! Do Hoschek parameter correction
  use precision
  implicit none
  ! Input/Output
  real(kind=realType)  ,intent(in)      :: tu(ku+nctlu),tv(kv+nctlv)
  real(kind=realType)  ,intent(inout)   :: u(nv,nu),v(nv,nu)
  real(kind=realType)  ,intent(in)      :: coef(ndim,nctlv,nctlu)
  integer           ,intent(in)      :: ku,kv,nctlu,nctlv,ndim,nu,nv
  real(kind=realType)  ,intent(in)      :: X(ndim,nv,nu)
 
  ! Working
  integer                            :: i,j,idim
  real(kind=realType)                   :: val(ndim),D(ndim)
  real(kind=realType)                   :: compute_rms_surface
  compute_rms_surface = 0.0
  do i=1,nu
     do j=1,nv
        call eval_surface(u(j,i),v(j,i),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
        D = val-X(:,j,i)
        do idim=1,ndim
           compute_rms_surface = compute_rms_surface + D(idim)**2
        end do
     end do
  end do
  compute_rms_surface = sqrt(compute_rms_surface/(nu*nv))

end function compute_rms_surface


! subroutine surface_jacobian_wrap2(u,v,tu,tv,ku,kv,nctlu,nctlv,nu,nv,Jac)

!   use precision
!   implicit none
!   ! Input
!   integer         , intent(in)          :: ku,kv,nctlu,nctlv,nu,nv
!   real(kind=realType), intent(in)          :: u(nu,nv),v(nu,nv)
!   real(kind=realType), intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
!   real(kind=realType), intent(out)         :: Jac(nu*nv,nctlu*nctlv)
!   ! Working
!   real(kind=realType)                      :: vniku(ku),worku(2*ku)
!   integer                               :: ilou,ileftu,mflagu

!   real(kind=realType)                      :: vnikv(kv),workv(2*kv)
!   integer                               :: ilov,ileftv,mflagv

!   integer                               :: i,j,ii,jj,iwork

!   ilou = 1
!   ilov = 1
!   do i=1,nu
!      do j = 1,nv
!         ! Get u interval
!         call intrv(tu,nctlu+ku,u(i,j),ilou,ileftu,mflagu)
!         if (mflagu == 0) then
!            call bspvn(tu,ku,ku,1,u(i,j),ileftu,vniku,worku,iwork)
!         else if (mflagu == 1) then
!            ileftu = nctlu
!            vniku(:) = 0.0
!            vniku(ku) = 1.0
!         end if

!         ! Get v interval
!         call intrv(tv,nctlv+kv,v(i,j),ilov,ileftv,mflagv)
!         if (mflagv == 0) then
!            call bspvn(tv,kv,kv,1,v(i,j),ileftv,vnikv,workv,iwork)
!         else if (mflagv == 1) then
!            ileftv = nctlv
!            vnikv(:) = 0.0
!            vnikv(kv) = 1.0
!         end if
        
!         do ii=1,ku
!            do jj = 1,kv
!               Jac( (i-1)*nv + j, (ileftu-ku+ii-1)*Nctlv + (ileftv-kv+jj)) = &
!                    vniku(ii)*vnikv(jj)
!            end do
!         end do
!      end do
!   end do
! end subroutine surface_jacobian_wrap2




! subroutine surface_para_corr(tu,tv,ku,kv,u,v,coef,nctlu,nctlv,ndim,nu,nv,X,rms)

!   ! Do Hoschek parameter correction
!   use precsion
!   implicit none
!   ! Input/Output
!   real(kind=realType)  ,intent(in)      :: tu(ku+nctlu),tv(kv+nctlv)
!   real(kind=realType)  ,intent(inout)   :: u(nu,nv),v(nu,nv)
!   real(kind=realType)  ,intent(in)      :: coef(nctlu,nctlv,ndim)
!   integer           ,intent(in)      :: ku,kv,nctlu,nctlv,ndim,nu,nv
!   real(kind=realType)  ,intent(in)      :: X(nu,nv,ndim)
!   real(kind=realType)  ,intent(out)     :: rms
!   ! Working
!   integer                            :: i,j,jj,max_inner_iter
!   real(kind=realType)                   :: lengthu,lengthv
!   real(kind=realType)                   :: D(ndim),D2(ndim),Dnorm,D2norm
!   real(kind=realType)                   :: val(ndim),deriv(2,ndim)
!   real(kind=realType)                   :: delta_c,delta_d,u_tilde,v_tilde
!   integer                            :: adj_u,adj_v
!   !Functions
!   real(kind=realType)                   :: norm,poly_length

!   max_inner_iter = 10
!   rms = 0.0
!   do i=1,nu
!      do j = 1,nv

!         lengthu = poly_length(X(:,j,:),nu,ndim)
!         lengthv = poly_length(X(i,:,:),nv,ndim)

!         adj_u = 1 ! Adjust them by default
!         adj_v = 1

!         if (i==1 .or. i==nu) then
!            adj_u = 0
!         end if

!         if (j==1 .or. j == nv) then
!            adj_v = 0
!         end if

!         call eval_surface(u(i,j),v(i,j),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
!         call eval_surface_deriv(u(i,j),v(i,j),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv)

!         D = X(i,j,:)-val

!         if (adj_u == 1) then
!            delta_c = dot_product(D,deriv(1,:))/norm(deriv(1,:),ndim)

!            inner_loop1: do jj=1,max_inner_iter
!               u_tilde = u(i,j)+ delta_c*(tu(nctlu+ku)-tu(1))/lengthu
!               call eval_surface(u_tilde,v(i,j),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
!               D2 = X(i,j,:)-val
!               if (norm(D,ndim) .ge. norm(D2,ndim)) then
!                  u(i,j) = u_tilde
!                  exit inner_loop1
!               else
!                  delta_c = delta_c*0.5
!               end if
!            end do inner_loop1
!         end if

!         if (adj_v == 1) then
!            delta_d = dot_product(D,deriv(2,:))/norm(deriv(2,:),ndim)

!            inner_loop2: do jj=1,max_inner_iter
!               v_tilde = v(i,j)+ delta_d*(tv(nctlv+kv)-tv(1))/lengthv
!               call eval_surface(u(i,j),v_tilde,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
!               D2 = X(i,j,:)-val
!               if (norm(D,ndim) .ge. norm(D2,ndim)) then
!                  v(i,j) = v_tilde
!                  exit inner_loop2
!               else
!                  delta_d = delta_d*0.5
!               end if
!            end do inner_loop2
!         end if

!      end do
!   end do

!   ! Lets redo the full RMS
!   rms = 0.0
!   do i=1,nu
!      do j=1,nv
!         call eval_surface(u(i,j),v(i,j),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
!         D = X(i,j,:)-val
!         rms = rms + dot_product(D,D)
!      end do
!   end do
!   rms = sqrt(rms/(nu*nv))
  

! end subroutine surface_para_corr

