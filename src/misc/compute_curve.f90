subroutine compute_curve(s,X,t,k,n,nctl,ndim,coef)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: compute_curve is the main function for the
  !               generation fo NURBS curves It does both
  !               interpolating and LMS fits, variable NURBS weights
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
  !
  !     Ouput coef - Real,Array of NURBS coefficients and
  !     weights. Size (nctl,ndim+1)

  implicit none
  
  ! Input
  integer         , intent(in)          :: k,nctl,ndim,n
  double precision, intent(in)          :: X(n,ndim)
  double precision, intent(inout)       :: s(n)
  double precision, intent(in)          :: t(nctl+k)

  ! Output
  double precision, intent(out)         :: coef(nctl,ndim+1)

  ! Working
  integer                               :: i,idim,inbv,iter,nrow
  integer                               :: niter,rank
  double precision                      :: work(3*k)
  double precision                      :: weight
  double precision                      :: Xcopy(n,ndim)
  integer                               :: info,lwork,ncol
  double precision                      :: length,res,norm,rms,val(ndim)

  double precision ,allocatable,dimension(:) :: work2,Svd
  double precision ,allocatable,dimension(:,:) :: Jac
  double precision ,allocatable,dimension(:) :: rhs
  
  ! Functions called
  double precision                      :: compute_lms_norm,poly_length
  
  print *,'Welcome to Compute Curve'
  !Initialization
  niter = 1
  coef(:,:) =  0.0
  coef(:,ndim+1) = 1.0 ! Initialize all the weights to 1
  length = poly_length(X,n,ndim)  

  ! Determine the work size (for the non-linear case)
  lwork = -1
  ncol = nctl*(ndim+1)-1
  nrow = n*ndim
  call DGELSS(nrow,ncol,ndim,Jac,nrow,Xcopy,max(nrow,ncol),Svd,-1,rank,work,lwork,info)
  lwork = work(1)
  allocate(work2(lwork))
  print *,'lwork',lwork
  ! Do one LMS fit to get the Initial Starting Point
  ncol = nctl
  nrow = n
  allocate(Jac(nrow,ncol),Svd(min(nrow,ncol)))
  call curve_jacobian_linear(Jac,nrow,ncol,ndim,t,k,s)
  Xcopy(:,:) = X(:,:)
  call DGELSS(nrow,ncol,ndim,Jac,nrow,Xcopy,max(nrow,ncol),Svd,-1,rank,work2,lwork,info)
  coef(:,1:2) = Xcopy(1:ncol,:)
  deallocate(Jac,Svd)

  rms = sqrt(compute_lms_norm(Xcopy,nrow,ncol,ndim)/n)  
  print *,'Linear LMS RMS:',rms
  
  ! Now do non-linear LMS fitting for NURBS with parameter correction
  ncol = nctl*(ndim+1)
  nrow = n*ndim
  
  allocate(Jac(nrow,ncol),Svd(min(nrow,ncol)),rhs(nrow))
  call curve_jacobian_non_linear(Jac,nrow,ncol,n,nctl,ndim,t,k,s,coef)
  print *,'nrow,ncol:',nrow,ncol
  !print *,'Jac',Jac
   do iter=1,niter
      print *,'iter:',iter
     ! Compute the RHS
      do i=1,n
         call eval_curve(s(i),t,k,coef,nctl,ndim,val  )
         do idim=1,ndim
            rhs((i-1)*ndim + idim) = X(i,idim)-val(idim)
         end do
      end do
      print *,'rhs:',rhs
     
      call DGELSS(nrow,ncol,1,Jac,nrow,rhs,max(nrow,ncol),Svd,-1,rank,work2,lwork,info)
      !    DGELSS(M,  N, NRHS,A  ,LDA,B  ,LDB, S,RCOND,RANK,WORK,LWORK,INFO)
      ! Update the coefficients

      do i=1,nctl
         do idim=1,ndim+1
            coef(i,idim) = coef(i,idim) + rhs((i-1)*(ndim+1)+idim)
            print *,'i,idim,delta:',i,idim,rhs((i-1)*(ndim+1)+idim)
         end do
      end do

      coef(:,ndim+1) = 1.0
     rms = sqrt(compute_lms_norm(rhs,nrow,ncol,1)/(nrow))
     print *,'Non-Linear LMS RMS:',rms
     
     !call curve_para_corr(t,k,coef,nctl,ndim,s,length,n,X)
     !call curve_jacobian_non_linear(Jac,nrow,ncol,n,nctl,ndim,t,k,s,coef)
   end do
  print *,'done all'
     
  deallocate(work2)
  deallocate(Jac,Svd,rhs)
end subroutine compute_curve

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


subroutine curve_jacobian_linear(Jac,n,nctl,ndim,t,k,s)

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

end subroutine curve_jacobian_linear

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


subroutine curve_para_corr(t,k,coef,nctl,ndim,s,length,n,X)
  ! Do Hoschek parameter correction
  implicit none
  double precision                   :: t(k+nctl)
  double precision                   :: coef(nctl,ndim+1)
  double precision                   :: s(n)
  double precision                   :: length
  double precision                   :: X(n,ndim)
  integer                            :: k,nctl,ndim,n,i,j,max_inner_iter

  double precision                   :: D(ndim),D2(ndim),val(ndim),deriv(ndim)
  double precision                   :: c,s_tilde,norm

  max_inner_iter = 10

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
           exit inner_loop
        else
           c = c*0.5
        end if
     end do inner_loop
 
  end do
end subroutine curve_para_corr

! subroutine compute_curve(s,X,t,k,n,nctl,ndim,coef,Jac)

!   !***DESCRIPTION
!   !
!   !     Written by Gaetan Kenway
!   !
!   !     Abstract: compute_curve is the main function for the
!   !               generation fo NURBS curves It does both
!   !               interpolating and LMS fits, variable NURBS weights
!   !               as well as Hoschek's Parameter Correction.
!   !
!   !     Description of Arguments
!   !     Input
!   !     s       - Real, Vector of s coordinates, length n
!   !     X       - Real, Array of X values to fit, Size (n,ndim)
!   !     t       - Real,Knot vector. Length nctl+k
!   !     k       - Integer,order of spline
!   !     nctl    - Integer,Number of control points
!   !     ndim    - Integer, spatial dimension of curve
!   !
!   !     Ouput coef - Real,Array of NURBS coefficients and
!   !     weights. Size (nctl,ndim+1)

!   implicit none
  
!   ! Input
!   integer         , intent(in)          :: k,nctl,ndim,n
!   double precision, intent(in)          :: X(n,ndim)
!   double precision, intent(inout)       :: s(n)
!   double precision, intent(in)          :: t(nctl+k)

!   ! Output
!   double precision, intent(out)         :: coef(nctl,ndim+1)

!   ! Working
!   integer                               :: i,idim,inbv,iter
!   integer                               :: niter,rank
!   double precision                      :: work(3*k)
!   double precision                      :: weight
!   double precision                      :: Xcopy(n,ndim)
!   double precision   ,intent(out)       :: Jac(n,nctl)
!   double precision                      :: Svd(nctl)
!   integer                               :: info,lwork
!   double precision                      :: length,res,norm,tot,val(ndim)

!   double precision ,allocatable,dimension(:) :: work2
  
!   ! Functions called
!   double precision                      :: compute_lms_norm,poly_length

!   print *,'Welcome to Compute Curve'

!   ! Determine the work size
!   lwork = -1
!   call DGELSS(n,nctl,ndim,Jac,n,Xcopy,n, Svd, -1 , rank,work,lwork,info)
!   lwork = work(1)
!   allocate(work2(lwork))
     
!   !Initialization
!   niter = 1000
!   coef(:,:) =  0.0
!   coef(:,ndim+1) = 1.0 ! Initialize all the weights to 1

!   length = poly_length(X,n,ndim)  

!   call curve_jacobian(Jac,n,nctl,ndim,t,k,s,coef)
  
!   do iter=1,niter
!      Xcopy(:,:) = X(:,:)
!      call DGELSS(n,nctl,ndim,Jac,n,Xcopy,n, Svd, -1 , rank,work2,lwork,info)
!      !    DGELSS(M,  N, NRHS,A  ,LDA,B  ,LDB, S,RCOND,RANK,WORK,LWORK,INFO)

!      do i=1,nctl
!         do idim=1,ndim
!            coef(i,idim) = Xcopy(i,idim)
!         end do
!      end do
     
!      call curve_para_corr(t,k,coef,nctl,ndim,s,length,n,X)
     
!      ! Check the actual deltas
!      tot = 0.0
!      do i=1,n
!         call eval_curve(s(i),t,k,coef,nctl,ndim,val)
!         tot = tot + (val(1)-X(i,1))**2  + (val(2)-X(i,2))**2
!      end do

!      if (mod(iter,100) == 0) then
!         print *,'iter,rms:',iter,sqrt(tot/n)
!      end if
!      call curve_jacobian(Jac,n,nctl,ndim,t,k,s,coef)     
!   end do
!   deallocate(work2)
! end subroutine compute_curve


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
