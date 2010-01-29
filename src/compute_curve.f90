subroutine compute_curve(s,X,t,k,n,nctl,ndim,coef,niter,tol)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: compute_curve is the main function for the
  !               generation fo B-spline curves It does both
  !               interpolating and LMS fits, variable B-spline weights
  !               as well as Hoschek's Parameter Correction.
  !
  !     Description of Arguments
  !     Input
  !     s       - Real, Vector of s coordinates, length n
  !     X       - Real, Array of X values to fit, Size (n,ndim)
  !     t       - Real,Knot vector. Length nctl+k
  !     k       - Integer,order of spline
  !     nctl    - Integer,Number of control points
  !     ndim    - Integer, spatial dimension of curve
  !     niter   - Integer, number of hoscek para-correction iterations
  !     tol     - Integer, relative tolerance for convergence
  !
  !     Ouput coef - Real,Array of B-spline coefficients. Size (nctl,ndim)

  use lms_jacobian
  use  lsqrModule,        only : LSQR
  use  lsqrCheckModule,   only : Acheck, xcheck

  implicit none
  
  ! Input
  integer         , intent(in)          :: k,nctl,ndim,n
  double precision, intent(in)          :: X(n,ndim)
  double precision, intent(inout)       :: s(n)
  double precision, intent(inout)       :: t(nctl+k)
  integer         , intent(in)          :: niter
  double precision, intent(in)          :: tol

  ! Output
  double precision, intent(out)         :: coef(nctl,ndim)

  ! Working
  integer                               :: i,idim,iter,nout
  double precision                      :: weight
  double precision                      :: length,res,norm,tot,val(ndim)
  integer                               :: istop,itn
  double precision                      :: Anorm,Acond,rnorm, Arnorm,xnorm
  ! Functions called
  double precision                      :: poly_length,floor

  print *,'Welcome to Compute Curve'
     
  !Initialization
  coef(:,:) =  0.0
  length = poly_length(X,n,ndim)  

  ! Setup the jacobain module
  call setup_jacobian(n,nctl,k)
  
  call curve_jacobian_linear(t,k,s,n,nctl)

  do iter=1,niter
     ! Solve
     idim = 1
     !do idim=1,ndim
     print *,'calling lsqr'
     print *,'n:',n
     print *,'Nctl:',Nctl
     !print *,'X(:,idim)',X(:,idim)
     !print *,'col_ind:',col_ind
     !print *,'A_ptr:',row_ptr
     !print *,'vals:',vals
     
     nout   = 6
     open(nout,file='LSQR.txt',status='unknown')
     call Acheck(n,Nctl,Aprod1,Aprod2,nout,istop)
     print *,'istop',istop


      call LSQR( n, Nctl, Aprod1, Aprod2,X(:,idim),0.0, .False., &
           coef(:,idim), vals, 1e-12, 1e-12, 1e8, 1000,nout,&
           istop, itn, Anorm, Acond, rnorm, Arnorm, xnorm )
   !end do
      close(nout)
     print *,'Results:'
     print *,'istop:',istop
     print *,'itn:',itn
     print *,'Anorm:',Anorm
     print *,'Acond:',Acond
     print *,'rnorm:',rnorm
     print *,'Arnorm:',Arnorm
     print *,'xnorm:',xnorm

     !call curve_para_corr(t,k,s,coef,nctl,ndim,length,n,X,rms)
     !call curve_jacobian_linear
     

     print *,'Done LSQR'
     ! Do convergence Check to break early
     
  end do
  call kill_jacobian()

end subroutine compute_curve


subroutine curve_jacobian_linear(t,k,s,n,nctl)
  
  use lms_jacobian

  implicit none
  ! Input
  double precision     ,intent(in)      :: t(nctl+k)
  double precision     ,intent(in)      :: s(n)
  integer              ,intent(in)      :: k,nctl,n
  ! Output
  !double precision     ,intent(out)     :: vals(n*k)
  !integer              ,intent(out)     :: col_ind(n*k),row_ptr(n)
  ! Working
  double precision                      :: vnikx(k),work(2*k)
  integer                               :: i,j,counter
  integer                               :: ilo,ileft,mflag,iwork
  
  ilo = 1
  counter = 1
  do i=1,n
     call intrv(t,nctl+k,s(i),ilo,ileft,mflag)
     if (mflag == 0) then
        call bspvn(t,k,k,1,s(i),ileft,vnikx,work,iwork)
     else if (mflag == 1) then
        ileft = nctl
        vnikx(:) = 0.0
        vnikx(k) = 1.0
     end if

     row_ptr(i) = counter
     do j=1,k
        col_ind(counter) = ileft-k+j
        vals(counter) =vnikx(j)
        counter = counter + 1
     end do

  end do
  row_ptr(n+1) = counter
end subroutine curve_jacobian_linear

function poly_length(X,n,ndim)
  ! Compute the length of the spatial polygon
  implicit none
  !Input
  double precision ,intent(in)    :: X(n,ndim)
  integer          ,intent(in)    :: n,ndim
  integer                         :: i,idim
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



subroutine curve_para_corr(t,k,s,coef,nctl,ndim,length,n,X,rms)

  ! Do Hoschek parameter correction
  implicit none
  ! Input/Output
  double precision  ,intent(in)      :: t(k+nctl)
  double precision  ,intent(inout)   :: s(n)
  double precision  ,intent(in)      :: coef(nctl,ndim)
  integer           ,intent(in)      :: k,nctl,ndim,n
  double precision  ,intent(in)      :: X(n,ndim)
  double precision  ,intent(in)      :: length
  double precision  ,intent(out)     :: rms
  ! Working
  integer                            :: i,j,max_inner_iter
  double precision                   :: D(ndim),D2(ndim),Dnorm,D2norm
  double precision                   :: val(ndim),deriv(ndim)
  double precision                   :: c,s_tilde,norm

  max_inner_iter = 10
  rms = 0.0
  do i=2,n-1
     call eval_curve(s(i),t,k,coef,nctl,ndim,val)
     call eval_curve_deriv(s(i),t,k,coef,nctl,ndim,deriv)
          
     D = X(i,:)-val
     c = dot_product(D,deriv)/norm(deriv,ndim)

     inner_loop: do j=1,max_inner_iter

        s_tilde = s(i)+ c*(t(nctl+k)-t(1))/length
        call eval_curve(s_tilde,t,k,coef,nctl,ndim,val)
        D2 = X(i,:)-val
        if (norm(D,ndim) .ge. norm(D2,ndim)) then
           s(i) = s_tilde
           rms = rms + dot_product(D2,D2)
           exit inner_loop
        else
           c = c*0.5
           if (j==max_inner_iter) then
              rms = rms + dot_product(D2,D2)
           end if
        end if
     end do inner_loop
  end do

  ! Add the first and Last
  call eval_curve(s(1),t,k,coef,nctl,ndim,val)
  rms = rms + dot_product(X(1,:)-val,X(1,:)-val)
  call eval_curve(s(n),t,k,coef,nctl,ndim,val)
  rms = rms + dot_product(X(n,:)-val,X(n,:)-val)
   rms = sqrt(rms/n)
end subroutine curve_para_corr



! subroutine curve_jacobian(Jac,n,nctl,ndim,t,k,s,coef)

!   implicit none
!   integer                               :: k,nctl,ndim,n
!   double precision                      :: Jac(n,nctl)
!   double precision                      :: s(n)
!   double precision                      :: t(nctl+k)
!   double precision                      :: coef(nctl,ndim+1)
!   double precision                      :: vnikx(k),work(2*k)
!   integer                               :: i,j,idim
!   integer                               :: ilo,ileft,mflag,iwork
  
!   Jac(:,:) = 0.0
!   ilo = 1
!   do i=1,n
!      call intrv(t,nctl+k,s(i),ilo,ileft,mflag)
!      if (mflag == 0) then
!         call bspvn(t,k,k,1,s(i),ileft,vnikx,work,iwork)
!      else if (mflag == 1) then
!         ileft = nctl
!         vnikx(:) = 0.0
!         vnikx(k) = 1.0
!      end if

!      do j=1,k
!         Jac(i,ileft-k+j) = vnikx(j)
!      end do
!   end do

! end subroutine curve_jacobian


function compute_lms_norm(X,m,n,ndim)
  ! Compute the norm from the dgels calc in general
  ! m is rows
  ! n is cols (m>n)
  ! ndim is spatial dimension
  implicit none
  integer                         :: m,n,ndim,i,idim
  double precision                :: X(m,ndim)
  double precision                :: compute_lms_norm

  compute_lms_norm = 0.0
  do i=n+1,m
     do idim=1,ndim
        compute_lms_norm = compute_lms_norm + X(i,idim)**2
     end do
  end do
end function compute_lms_norm



subroutine curve_jacobian_linear2(Jac,n,nctl,ndim,t,k,s)

  implicit none
  integer                               :: k,nctl,ndim,n
  double precision                      :: Jac(n,nctl)
  double precision                      :: s(n)
  double precision                      :: t(nctl+k)
  double precision                      :: vnikx(k),work(2*k)
  integer                               :: i,j
  integer                               :: ilo,ileft,mflag,iwork
  
  Jac(:,:) = 0.0
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

end subroutine curve_jacobian_linear2

function norm(X,n)
  ! Compute the L2 nomr of X
  implicit none
  double precision       :: X(n)
  double precision       :: norm
  integer                :: i,n
  norm = 0.0
  do i=1,n
     norm = norm + X(i)**2
  end do
  norm = sqrt(norm)
end function norm

subroutine curve_jacobian_non_linear(Jac,nrow,ncol,n,nctl,ndim,t,k,s,coef)

  implicit none
  integer                               :: k,nrow,ncol,n,nctl,ndim
  double precision                      :: Jac(nrow,ncol)
  double precision                      :: t(nctl+k)
  double precision                      :: s(n)
  double precision                      :: coef(nctl,ndim+1)
  double precision                      :: val(n,ndim),val0(n,ndim)
  integer                               :: idim,jdim,i,j
  Jac(:,:) = 0.0
  ! FD the whole fucking thing

  call eval_curve_V(s,t,k,coef,nctl,ndim,n,val0)
  do j=1,nctl               !-> Loop over the Cols
     do idim=1,ndim+1       !- > Loop over the COls
        
        coef(j,idim) = coef(j,idim) + 1e-5
        call eval_curve_V(s,t,k,coef,nctl,ndim,n,val)
        
        ! Set the stuff
        do i = 1,n
           do jdim=1,ndim
              Jac(ndim*(i-1) + jdim  ,(j-1)*(ndim+1) + idim)  = &
                   (val(i,jdim) -val0(i,jdim))/1e-6
           end do
        end do
        
        coef(j,idim) = coef(j,idim) - 1e-5

     end do
  end do

end subroutine curve_jacobian_non_linear




! subroutine curve_jacobian_non_linear(Jac,nrow,ncol,n,nctl,ndim,t,k,s,coef)

!   implicit none
!   integer                               :: k,nrow,ncol,n,nctl,ndim
!   double precision                      :: Jac(nrow,ncol)
!   double precision                      :: t(nctl+k)
!   double precision                      :: s(n)
!   double precision                      :: coef(nctl,ndim+1)
!   double precision                      :: vnikx(k),work(3*k)
!   double precision                      :: val(ndim),deriv(ndim),val2(ndim)
!   double precision                      :: R,weight
!   integer                               :: i,j,idim,icol
!   integer                               :: inbv
!   integer                               :: ilo,ileft,mflag,iwork
  
!   double precision                      :: bvalu
!   Jac(:,:) = 0.0
!   ilo = 1
!   inbv = 1
!   do i=1,n

!      call eval_curve      (s(i),t,k,coef,nctl,ndim,val  )
!      call eval_curve_deriv(s(i),t,k,coef,nctl,ndim,deriv)
!      weight = bvalu(t,coef(:,ndim+1),nctl,k,0,s(i),inbv,work)
     
!      call intrv(t,nctl+k,s(i),ilo,ileft,mflag)
!      if (mflag == 0) then
!         call bspvn(t,k,k,1,s(i),ileft,vnikx,work,iwork)
!      else if (mflag == 1) then
!         ileft = nctl
!         vnikx(:) = 0.0
!         vnikx(k) = 1.0
!      end if

!      ! Each Block should look like this: x,y,z for a point row wise
!      ! and Cx,Cy,Cz,weight for a single control point along the columns
!      ! [ x  0  0 x ]
!      ! [ 0  x  0 x ]
!      ! [ 0  0  x x ]
!      do j=1,k
!         icol = (ileft-k+j) ! Control Point of interest
!         ! Do the Rational Basis First Entries First
!         R = vnikx(j)*coef(icol,ndim+1)/weight

!         if (icol > 1 .and. icol<ncol-1) then
!            do idim=1,ndim
!               Jac(ndim*(i-1) + idim, (ndim+1)*(icol-2) +ndim + idim) = R
!               print *,'i,j,val:',ndim*(i-1)+idim,(ndim+1)*(icol-2)+ndim+idim,R
!            end do
!         else if (icol == 1) then
!            do idim=1,ndim
!               Jac(ndim*(i-1) + idim, idim) = R
!               print *,'i,j,val:',ndim*(i-1)+idim,idim,R
!            end do
!         end if

!        !  ! Next the weight derivatives with fd
!          coef(icol,ndim+1) = coef(icol,ndim+1) + 1e-6

!          call eval_curve(s(i),t,k,coef,nctl,ndim,val2)
!          deriv = (val2-val)/1e-6
!          coef(icol,ndim+1) = coef(icol,ndim+1) - 1e-6

! !         print *,'fd'
         
!          if (icol>1 .and. icol<nctl-1) then
!             do idim=1,ndim
!                Jac(ndim*(i-1) + idim,(ndim+1)*(icol-2) + ndim + ndim+1) = deriv(idim)
!                print *,'i,j,val:',ndim*(i-1)+idim,(ndim+1)*(icol-1)+ndim+1,deriv(idim)
!             end do
!          end if
!           !print *,'analytic'
!           !do idim =1,ndim
!            !Jac(ndim*(i-1) + idim,(ndim+1)*(icol-1) + ndim+1) = &
!             !    vnikx(j)*(coef(icol,idim) - val(idim))/weight
!            !print *,'i,j,val:',ndim*(i-1)+idim,(ndim+1)*(icol-1)+ndim+1,vnikx(j)*(coef(icol,idim) - val(idim))/weight
!        !end do
!        end do


!   end do

! end subroutine curve_jacobian_non_linear
