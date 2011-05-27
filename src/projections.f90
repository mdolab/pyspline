subroutine point_curve(x0,t,k,coef,nctl,ndim,N,Niter,eps1,eps2,s,Diff)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: point_curve attempts to solve the point inversion problem
  !
  !     Description of Arguments
  !     Input
  !     x0      - Real, array size(ndim,N) N points we are trying to invert
  !     t       - Real, Knot vector. Length nctl+k
  !     k       - Integer,order of B-spline
  !     coef    - Real,Array of B-spline coefficients and weights. Size (ndim,nctl)
  !     nctl    - Integer,Number of control points
  !     ndim    - Integer, spatial dimension of curve
  !     Niter   - Integer, Maximum number of Netwton iterations
  !     eps1    - Real - Eculdian Distance Convergence Measure
  !     eps2    - Real - Cosine Convergence Measure
  !     s       - Real,vector, length(N), guess parameters where C(s)
  !               is closest to x0 
  !
  !     Ouput 
  !     s       - Real, vector, length(N) parameter where C(s) is closest to x0
  !     diff    - Real,array, size(ndim,N)- Distance between x0 and curve(s)

  implicit none
  ! Input
  double precision, intent(in)          :: x0(ndim,N)
  integer         , intent(in)          :: k,nctl,ndim,N,niter
  double precision, intent(in)          :: t(nctl+k)
  double precision, intent(in)          :: coef(ndim,nctl)
  double precision, intent(in)          :: eps1,eps2

  ! Output
  double precision, intent(out)         :: s(N),Diff(ndim,N)

  ! Working
  double precision                      :: val(ndim),deriv(ndim),deriv2(ndim)
  double precision                      :: val0(ndim),s0(N)
  integer                               :: i,j,ipt,counter,max_inner_iter
  integer                               :: istart,nss
  double precision                      :: D,D0,delta,D2(ndim)
  integer                               :: n_sub ! Huristic Value
  logical                               :: brute_force
  ! Alloctable
  double precision,allocatable          :: curve_vals(:,:),ss(:)
  double precision :: c1
  ! Functions
  double precision                      :: norm,ndp

  ! Initialization
  n_sub=1  ! Is 3  Good Here?? Huristic!
  brute_force = .False.
  max_inner_iter = 20
  ! Determine if any of the guesses are out of range, if so, do the brute force search

  point_loop: do ipt=1,N
     if (s(ipt) < 0 .or. s(ipt) > 1) then
        brute_force = .True.
        exit point_loop
     end if
  end do point_loop

  if (brute_force) then
     ! Dynamically allot the required space for the brute-force search search
         nss = nctl-k+2 + (nctl-k+1)*n_sub
     allocate (ss(nss))
     allocate(curve_vals(ndim,nss))

     counter = 1
     do i=1,nctl-k+1
        if (i==1) then
           istart = 0
        else
           istart = 1
        end if
        do j=istart,n_sub+1
           ss(counter) = t(i+k-1) + (real(j)/(n_sub+1)*(t(i+k)-t(i+k-1)))
           counter = counter + 1
        end do
         end do
  
     do i=1,nss
        call eval_curve(ss(i),t,k,coef,nctl,ndim,1,curve_vals(:,i))
     end  do

     s0(:) = s(:)
     ! Now do the actual searching
     do ipt=1,N
        
        if (s(ipt) < 0 .or. s(ipt) > 1) then ! Still only do it if we have a bad guess
           D0 = norm(x0(:,ipt)-curve_vals(:,1),ndim)
           s0(ipt) = ss(1)
           do i=1,nss
              D = norm(curve_vals(:,i)-x0(:,ipt),ndim)
              if (D<D0) then
                 s0(ipt) = ss(i)
                 D0 = D
              end if
           end do
        else
           s0(i) = s(i)
        end if
     end do
  else
     s0(:) = s(:)   
  end if

  s = s0
  do ipt=1,N
     call eval_curve(s0(ipt),t,k,coef,nctl,ndim,1,val)
     call eval_curve_deriv(s0(ipt),t,k,coef,nctl,ndim,deriv)
     call eval_curve_deriv2(s0(ipt),t,k,coef,nctl,ndim,deriv2)

     Diff(:,ipt) = val-x0(:,ipt)
    
     iteration_loop: do i=1,Niter
        ! Check the Convergence Criteria

        if (norm(Diff(:,ipt),ndim) <= eps1) then
               exit iteration_loop
        end if

        c1 = ndp(deriv,Diff(:,ipt),ndim)/(norm(deriv,ndim)*norm(Diff(:,ipt),ndim))
        if (c1 <= eps2) then
           exit iteration_loop
        end if

        s0(ipt) = s(ipt)
        call eval_curve(s0(ipt),t,k,coef,nctl,ndim,1,val)
        call eval_curve_deriv(s0(ipt),t,k,coef,nctl,ndim,deriv)
        call eval_curve_deriv2(s0(ipt),t,k,coef,nctl,ndim,deriv2)
        Diff(:,ipt) = val-x0(:,ipt)

        delta = -dot_product(deriv,Diff(:,ipt))/(dot_product(deriv2,Diff(:,ipt)) + norm(deriv,ndim)**2)

        ! Bounds Checking

        if (s0(ipt) + delta < t(1)) then
           delta = t(1)-s0(ipt)
        end if

        if (s0(ipt) + delta > t(nctl+k)) then
           delta = t(nctl+k) - s0(ipt)
        end if

        inner_loop: do j=1,max_inner_iter
          
           s(ipt) = s0(ipt) + delta
           
           call eval_curve(s(ipt),t,k,coef,nctl,ndim,1,val)
           D2 = val-x0(:,ipt)
           if ( ( norm(D2,ndim)-norm(diff(:,ipt),ndim)) .ge. eps1) then
              delta = delta * 0.5
           else
              exit inner_loop
           end if
        end do inner_loop

        if (norm((s(ipt)-s0(ipt))*deriv,ndim) <= eps1) then
           exit iteration_loop
        end if

     end do iteration_loop
  end do
  if (brute_force) then
     deallocate(curve_vals,ss)
  end if

end subroutine point_curve

subroutine point_surface(x0,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,N,niter,eps1,eps2,u,v,Diff)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: pint_surface attempts to solve the point inversion problem for a surface
  !
  !     Description of Arguments
  !     Input
  !     x0      - Real, array size(ndim,N) points we are trying to invert
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     coef    - Real,Array of B-spline coefficients  Size (ndim,nctlv,nctlu)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     ndim    - Integer, Spatial Dimension
  !     Niter   - Integer, Maximum number of Netwton iterations
  !     eps1    - Real - Eculdian Distance Convergence Measure
  !     eps2    - Real - Cosine Convergence Measure
  !
  !     Ouput 
  !     u       - Real,vector size(N), u parameters where S(u,v) is closest to x0
  !     v       - Real,vector size(N), v parameters where S(u,v) is closest to x0
  !     diff    - Real Array size(ndim,N) - Distance between x0 and S(u,v)

  implicit none
  ! Input
  double precision, intent(in)          :: x0(ndim,N)
  integer         , intent(in)          :: ku,kv,nctlu,nctlv,ndim,niter,N
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
  double precision, intent(in)          :: coef(ndim,nctlv,nctlu)
  double precision, intent(in)          :: eps1,eps2
  ! Output
  double precision, intent(inout)       :: u(N),v(N)
  double precision, intent(out)         :: diff(ndim,N)


  ! Working
  double precision                      :: val(ndim),deriv(ndim,2),deriv2(ndim,2,2)
  double precision                      :: val0(ndim)
  logical                               :: brute_force
  integer                               :: i,j,ii,jj,counter,ipt,max_inner_iter
  double precision                      :: D,D0,u0(N),v0(N),delta(2),D2(ndim)
  double precision                      :: A(2,2),ki(2)
  integer                               :: n_sub_u,n_sub_v,n_sub_w  ! Huristic Value
  integer                               :: istart,nuu,nvv

  ! Alloctable
  double precision,allocatable          :: surface_vals(:,:,:),uu(:),vv(:)

  ! Functions     
  double precision                      :: norm

  max_inner_iter = 20
  if (ku == 2) then
     n_sub_u = 2
  else
     n_sub_u = 2
  end if

  if (kv == 2) then
     n_sub_v = 2
  else
     n_sub_v = 2
  end if

  brute_force = .True.
  ! First we will evaluate the surface at n points inside each knot span in each direction

  !if we are given a bad guess do the brute force
  point_loop: do ipt=1,N
     if (u(ipt) < 0 .or. u(ipt) > 1 .or. v(ipt) < 0 .or. v(ipt) > 1) then
        brute_force = .True.
        exit point_loop
     end if
  end do point_loop
  if (brute_force) then
     ! Generate the uu and vv values
     nuu = nctlu-ku+2 + (nctlu-ku+1)*n_sub_u
     nvv = nctlv-kv+2 + (nctlv-kv+1)*n_sub_v
     allocate (uu(nuu),vv(nvv))
 
     counter = 1
     do i=1,nctlu-ku+1
        if (i==1) then
           istart = 0
        else
           istart = 1
        end if
        do j=istart,n_sub_u+1
           uu(counter) = tu(i+ku-1) + (real(j)/(n_sub_u+1)*(tu(i+ku)-tu(i+ku-1)))
           counter = counter + 1
        end do
     end do

     counter = 1
     do i=1,nctlv-kv+1
        if (i==1) then
           istart = 0
        else
           istart = 1
        end if
        do j=istart,n_sub_v+1
           vv(counter) = tv(i+kv-1) + (real(j)/(n_sub_v+1)*(tv(i+kv)-tv(i+kv-1)))
           counter = counter + 1
        end do
     end do
     allocate(surface_vals(ndim,nvv,nuu))

     do i=1,nuu
        do j=1,nvv
           call eval_surface(uu(i),vv(j),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,&
                1,1,surface_vals(:,j,i))
        end do
     end do

     ! Now do comparison
     u0(:) = u(:)   
     v0(:) = v(:)

     do ipt=1,N
        if (u(ipt) < 0 .or. u(ipt) > 1 .or. v(ipt) < 0 .or. v(ipt) > 1) then
           val0 = surface_vals(:,1,1)
           D0 = norm(val0-x0(ipt,:),ndim)
           u0(ipt) = uu(1)
           v0(ipt) = vv(1)
           do i = 1,nuu
              do j = 1,nvv
                 D = norm(surface_vals(:,j,i)-X0(:,ipt),ndim)
                 if (D<D0) then
                    u0(ipt) = uu(i)
                    v0(ipt) = vv(j)
                    D0 = D
                 end if
              end do
           end do
        else
           u0(ipt) = u(ipt)
           v0(ipt) = v(ipt)
        end if
     end do
  else
     u0(:) = u(:)   
     v0(:) = v(:)
  end if

  ! Now we have u0 and v0 so we can do the newton 
  u = u0
  v = v0

  do ipt=1,N

     call eval_surface(u(ipt),v(ipt),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,1,1,val)
     call eval_surface_deriv(u(ipt),v(ipt),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv)
     call eval_surface_deriv2(u(ipt),v(ipt),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv2)
     
     Diff(:,ipt) = val - x0(:,ipt)
     
     iteration_loop: do i=1,niter
        ! Check the convergence criteria
        if (norm(Diff(:,ipt),ndim) <= eps1) then
           exit iteration_loop
        end if

        if (norm(dot_product(deriv(:,1),Diff(:,ipt)),ndim)/(norm(deriv(:,1),ndim)*norm(Diff(:,ipt),ndim)) <= eps2 .and. &
             norm(dot_product(deriv(:,2),Diff(:,ipt)),ndim)/(norm(deriv(:,2),ndim)*norm(Diff(:,ipt),ndim)) <= eps2 ) then
           exit iteration_loop
        end if
        u0(ipt) = u(ipt)
        v0(ipt) = v(ipt)
        call eval_surface(u(ipt),v(ipt),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,1,1,val)
        call eval_surface_deriv(u(ipt),v(ipt),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv)
        call eval_surface_deriv2(u(ipt),v(ipt),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv2)

        Diff(:,ipt) = val-x0(:,ipt)

        A(1,1) = norm(deriv(:,1),ndim)**2 + dot_product(Diff(:,ipt),deriv2(:,1,1))
        A(1,2) = dot_product(deriv(:,1),deriv(:,2)) + dot_product(Diff(:,ipt),deriv2(:,1,2))
        A(2,1) = A(1,2)
        A(2,2) = norm(deriv(:,2),ndim)**2 + dot_product(Diff(:,ipt),deriv2(:,2,2))

        ki(1) = -dot_product(Diff(:,ipt),deriv(:,1))
        ki(2) = -dot_product(Diff(:,ipt),deriv(:,2))

        call solve_2by2(A,ki,delta)

        ! Bounds Checking
        if (u0(ipt)+delta(1) < tu(1)) then
           delta(1) = tu(1)-u0(ipt)
        end if

        if (u0(ipt)+delta(1) > tu(nctlu+ku)) then
           delta(1) = tu(nctlu+ku)-u0(ipt)
        end if

        if (v0(ipt)+delta(2) < tv(1)) then
           delta(2) = tv(1)-v0(ipt)
        end if

        if (v0(ipt)+delta(2) > tv(nctlv+kv)) then
           delta(2) = tv(nctlv+kv)-v0(ipt)
        end if

        inner_loop: do j=1,max_inner_iter
           u(ipt) = u0(ipt) + delta(1)
           v(ipt) = v0(ipt) + delta(2)
          
           call eval_surface(u(ipt),v(ipt),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,1,1,val)

           D2 = val-x0(:,ipt)
      
           if ( (norm(D2,ndim)-norm(diff(:,ipt),ndim)) .ge. eps1) then
              delta = delta * 0.5
           else
              exit inner_loop
           end if
        end do inner_loop
       
        ! No Change convergence Test

        if (norm( (u(ipt)-u0(ipt))*deriv(:,1) + (v(ipt)-v0(ipt))*deriv(:,2),ndim) <= eps1) then
           exit iteration_loop
        end if
     end do iteration_loop
  end do
  if (brute_force) then
     deallocate(surface_vals,uu,vv)
  end if
end subroutine point_surface

subroutine point_volume(x0,tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,N,&
     niter,eps,u,v,w,Diff)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: point_volume attempts to solve the point inversion problem for a volume 
  !
  !     Description of Arguments
  !     Input
  !     x0      - Real, array size(ndim,N) points we are trying to invert
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     tw      - Real,Knot vector in w. Length nctlw+kw
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     kw      - Integer, order of B-spline in w
  !     coef    - Real,Array of B-spline coefficients  Size (ndim,nctlw,nctlv,nctlu)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     nctlw   - Integer,Number of control points in w
  !     ndim    - Integer, Spatial Dimension
  !     Niter   - Integer, Maximum number of Netwton iterations
  !     eps     - Real - Eculdian Distance Convergence Measure
  !
  !     Ouput 
  !     u       - Real,vector size(N), u parameters where V(u,v,w) is closest to x0
  !     v       - Real,vector size(N), v parameters where V(u,v,w) is closest to x0
  !     w       - Real,vector size(N), w parameters where V(u,v,w) is closest to x0
  !     diff    - Real Array size(N,ndim) - Distance between x0 and V(u,v,w)

  implicit none
  ! Input
  double precision, intent(in)          :: x0(ndim,N)
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,ndim,niter,N
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  double precision, intent(in)          :: coef(ndim,nctlw,nctlv,nctlu)
  double precision, intent(in)          :: eps
  ! Output
  double precision, intent(out)         :: u(N),v(N),w(N)
  double precision, intent(out)         :: diff(ndim,N)

  ! Working
  double precision                      :: D,D0,u0(N),v0(N),w0(N),delta(3),D2(ndim)
  integer                               :: n_sub_u,n_sub_v,n_sub_w,nuu,nvv,nww
  integer                               :: istart,counter
  logical                               :: brute_force
  double precision                      :: val0(ndim)

  double precision                      :: val(ndim)
  double precision                      :: deriv(ndim,3)
  double precision                      :: deriv2(ndim,3,3)
  double precision                      :: R(nDim),low(nDim),high(nDim)
  double precision                      :: pt(nDim),newpt(nDim),update(nDim)
  double precision                      :: step,dist,nDist,pgrad
  double precision                      :: gnorm,grad_norm,wolfe
  double precision                      :: fval,nfval,c,p_diff


  double precision                      :: grad(nDim),hessian(nDim,nDim)
  integer                               :: iDim,jDim,ipt,i,j,k,m,nLine
  logical                               :: flag,cflag

  ! Alloctable
  double precision,allocatable          :: volume_vals(:,:,:,:)

  ! Functions     
  double precision                      :: norm

   ! Generate the uu and vv values
  nuu = nctlu-ku+2
  nvv = nctlv-kv+2 
  nww = nctlw-kw+2 

  allocate(volume_vals(ndim,nww,nvv,nuu))

  do i=1,nuu
     do j=1,nvv
        do k=1,nww
           call eval_volume(tu(ku-1+i),tv(kv-1+j),tw(kw-1+k),&
                tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,&
                1,1,1,volume_vals(:,k,j,i))
        end do
     end do
  end do
  ! Now do comparison
  u0(:) = u(:)   
  v0(:) = v(:)
  w0(:) = w(:)
  do ipt=1,N
     val0 = volume_vals(:,1,1,1)
     D0 = norm(val0-x0(ipt,:),ndim)
     u0(ipt) = tu(ku)
     v0(ipt) = tv(kv)
     w0(ipt) = tw(kw)
     do i = 1,nuu
        do j = 1,nvv
           do k = 1,nww
              
              D = norm(volume_vals(:,k,j,i)-X0(:,ipt),ndim)
              if (D<D0) then
                 u0(ipt) = tu(ku-1+i)
                 v0(ipt) = tv(kv-1+j)
                 w0(ipt) = tw(kw-1+k)
                 D0 = D
              end if
           end do
        end do
     end do
  end do

  ! Set lower and upper bounds for u,v,w based on knot vector
  low(1) = tu(1)
  low(2) = tv(1)
  low(3) = tw(1)
  high(1) = tu(Nctlu+ku)
  high(2) = tv(Nctlv+kv)
  high(3) = tw(Nctlw+kw)

  ! Number of line search iterations
  nLine = 20
  ! Tolerance for the strong wolfe line-search conditions
  wolfe = 0.001
  do ipt=1,N

     pt(1) = u0(ipt)
     pt(2) = v0(ipt)
     pt(3) = w0(ipt)

     iteration_loop: do i=1,niter
        call eval_volume(pt(1),pt(2),pt(3),tu,tv,tw,ku,kv,kw,coef,&
             nctlu,nctlv,nctlw,ndim,1,1,1,val)
        call eval_volume_deriv(pt(1),pt(2),pt(3),tu,tv,tw,ku,kv,kw,coef,&
             nctlu,nctlv,nctlw,ndim,deriv)
        call eval_volume_deriv2(pt(1),pt(2),pt(3),tu,tv,tw,ku,kv,kw,coef,&
             nctlu,nctlv,nctlw,ndim,deriv2)
        
        ! Distance is R, "function value" fval is what we minimize
        R = val - X0(:,ipt)
        nDist = norm(R,nDim)
        fval = 0.5*nDist**2
        
        ! Calculate the Gradient
        do idim=1,nDim
           grad(idim) = dot_product(R,deriv(:,idim))
        end do

        ! Calculate the Hessian
        do jDim=1,nDim
           do iDim=1,nDim
              hessian(iDim,jDim) = dot_product(deriv(:,iDim),deriv(:,jDim)) + &
                   dot_product(R,deriv2(:,iDim,jDim))
           end do
        end do

      !   ! Bounds Checking
        
!         do iDim=1,nDim
!            flag = .False.
!            if (pt(iDim) < low(iDim)+eps .and. grad(iDim) >= 0.0) then
!               flag = .True.
!               pt(iDim) = low(iDim)
!            end if

!            if (pt(iDim) > high(iDim)-eps .and. grad(iDim) <= 0.0) then
!               flag = .True.
!               pt(iDim) = high(iDim)
!            end if
!            if ( flag ) then
!               grad(iDim) = 0.0
!               hessian(:,iDim) = 0.0
!               hessian(iDim,:) = 0.0
!               hessian(iDim,iDim) = 1.0
!            end if
!         end do

        ! Check the norm of the gradient
        grad_norm = norm(grad,nDim)
        if (grad_norm < eps) then
           exit iteration_loop
        end if
        
        ! Invert the hessian, compute the update and the projected gradient
        call solve_3by3(hessian,grad,update)
        update = -update
        pgrad = dot_product(update,grad)

        !Check that this is a descent direction - 
        !otherwise use the negative gradient    
        if ( pgrad >= 0.0 ) then
           update = -grad/gnorm
           pgrad = dot_product(update,grad)
        end if

        step = 1.0
        nDist = 0.0
        lineloop: do m=1,nLine
          
           newpt = pt + step * update

           cflag = .False. ! Check if the constraints are applied
           ! Check if the new point exceeds the bounds
           do iDim=1,nDim
              if (newpt(iDim) > high(iDim)) then
                 cflag = .True.
                 newpt(iDim) = high(iDim)
              else if (newpt(iDim) < low(iDim)) then
                 cflag = .True.
                 newpt(iDim) = low(iDim)
              end if
           end do
           
           ! Re-evaluate new position
         
           call eval_volume(newpt(1),newpt(2),newpt(3),&
                tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,1,1,1,val)

           ! Check if the new function value is lower, 
           ! otherwise adjust the step size
           R = val - X0(:,ipt)
           ndist = norm(R,nDim)
           nfval = 0.5*ndist**2

           ! Check if the new point satisfies the wolfe condition
           if ( nfval < fval + pgrad * wolfe * step ) then
              dist = ndist
              exit lineloop
           end if

           ! Calculate the new step length
           if ( cflag ) then
              ! If the constraints are applied - and the new point
              ! doesn't satisfy the Wolfe conditions, it doesn't make
              ! sense to apply a quadratic approximation
              step = 0.25 * step
           else
              ! c = nfval - fval - pgrad * step is always positive since
              ! nfval - fval > pgrad * wolfe * step > pgrad * step
              c = ( ( nfval - fval ) - pgrad * step )
              step = - step * step * pgrad/( 2.0 * c ) 
              ! This update is always less than the original step length
           end if
        end do lineloop

        if ( m == nline ) then
           dist = ndist
        else
           ! Check if there has been no change in the coordinates
           p_diff = norm(pt-newpt,nDim)
           if (p_diff < eps) then
              exit Iteration_loop
           end if
        end if
        pt = newpt

     end do iteration_loop
     
     ! Set the final values of the parameters and our distance function
     u(ipt) = pt(1)  
     v(ipt) = pt(2)
     w(ipt) = pt(3)
     diff(:,ipt) = R
  end do

  deallocate(volume_vals)

end subroutine point_volume

! subroutine point_volume_old(x0,tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,N,niter,eps1,eps2,n_sub,u,v,w,Diff)

!   !***DESCRIPTION
!   !
!   !     Written by Gaetan Kenway
!   !
!   !     Abstract: point_surface attempts to solve the point inversion problem for a surface
!   !
!   !     Description of Arguments
!   !     Input
!   !     x0      - Real, array size(ndim,N) points we are trying to invert
!   !     tu      - Real,Knot vector in u. Length nctlu+ku
!   !     tv      - Real,Knot vector in v. Length nctlv+kv
!   !     tw      - Real,Knot vector in w. Length nctlw+kw
!   !     ku      - Integer, order of B-spline in u
!   !     kv      - Integer, order of B-spline in v
!   !     kw      - Integer, order of B-spline in w
!   !     coef    - Real,Array of B-spline coefficients  Size (ndim,nctlw,nctlv,nctlu)
!   !     nctlu   - Integer,Number of control points in u
!   !     nctlv   - Integer,Number of control points in v
!   !     nctlw   - Integer,Number of control points in w
!   !     ndim    - Integer, Spatial Dimension
!   !     Niter   - Integer, Maximum number of Netwton iterations
!   !     eps1    - Real - Eculdian Distance Convergence Measure
!   !     eps2    - Real - Cosine Convergence Measure
!   !     n_sub   - Integer - Overide for projections
!   !
!   !     Ouput 
!   !     u       - Real,vector size(N), u parameters where V(u,v,w) is closest to x0
!   !     v       - Real,vector size(N), v parameters where V(u,v,w) is closest to x0
!   !     w       - Real,vector size(N), w parameters where V(u,v,w) is closest to x0
!   !     diff    - Real Array size(N,ndim) - Distance between x0 and V(u,v,w)

!   implicit none
!   ! Input
!   double precision, intent(in)          :: x0(ndim,N)
!   integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,ndim,niter,N
!   double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
!   double precision, intent(in)          :: coef(ndim,nctlw,nctlv,nctlu)
!   double precision, intent(in)          :: eps1,eps2
!   ! Output
!   double precision, intent(inout)       :: u(N),v(N),w(N)
!   double precision, intent(out)         :: diff(ndim,N)


!   ! Working
!   double precision                      :: val(ndim)
!   double precision                      :: deriv(ndim,3)
!   double precision                      :: deriv2(ndim,3,3)
!   double precision                      :: val0(ndim)
!   logical                               :: brute_force
!   integer                               :: i,j,k,ii,jj,kk,counter,ipt,max_inner_iter
!   double precision                      :: D,D0,u0(N),v0(N),w0(N),delta(3),D2(ndim)
!   double precision                      :: A(3,3),ki(3)
!   integer                               :: n_sub_u,n_sub_v,n_sub_w,n_sub ! Huristic Value
!   integer                               :: istart,nuu,nvv,nww

!   ! Alloctable
!   double precision,allocatable          :: volume_vals(:,:,:,:)
!   double precision,allocatable          :: uu(:),vv(:),ww(:)
!   double precision :: c1,c2,c3
!   ! Functions     
!   double precision                      :: norm,dotproduct,ndp

!   max_inner_iter = 20

!   if (n_sub < 0) then
!      if (ku == 2) then
!         n_sub_u = 5
!      else
!         n_sub_u = 2
!      end if

!      if (kv == 2) then
!         n_sub_v = 5
!      else
!         n_sub_v = 2
!      end if

!      if (kw == 2) then
!         n_sub_w = 5
!      else
!         n_sub_w = 2
!      end if
!   else
!      n_sub_u = 5!n_sub
!      n_sub_v = 5!n_sub
!      n_sub_w = 5!n_sub
!   end if

!   brute_force = .True.

!   !if we are given a bad guess do the brute force
!   point_loop: do ipt=1,N
!      if (u(ipt) < 0 .or. u(ipt) > 1 .or. v(ipt) < 0 .or. v(ipt) > 1) then
!         brute_force = .True.
!         exit point_loop
!      end if
!   end do point_loop

!   if (brute_force) then
!    ! Generate the uu and vv values
!      nuu = nctlu-ku+2 + (nctlu-ku+1)*n_sub_u
!      nvv = nctlv-kv+2 + (nctlv-kv+1)*n_sub_v
!      nww = nctlw-kw+2 + (nctlw-kw+1)*n_sub_w
!      allocate (uu(nuu),vv(nvv),ww(nww))
     

!      counter = 1
!      do i=1,nctlu-ku+1
!         if (i==1) then
!            istart = 0
!         else
!            istart = 1
!         end if
!         do j=istart,n_sub_u+1
!            uu(counter) = tu(i+ku-1) + (real(j)/(n_sub_u+1)*(tu(i+ku)-tu(i+ku-1)))
!            counter = counter + 1
!         end do
!      end do

!      counter = 1
!      do i=1,nctlv-kv+1
!         if (i==1) then
!            istart = 0
!         else
!            istart = 1
!         end if
!         do j=istart,n_sub_v+1
!            vv(counter) = tv(i+kv-1) + (real(j)/(n_sub_v+1)*(tv(i+kv)-tv(i+kv-1)))
!            counter = counter + 1
!         end do
!      end do

!      counter = 1
!      do i=1,nctlw-kw+1
!         if (i==1) then
!            istart = 0
!         else
!            istart = 1
!         end if
!         do j=istart,n_sub_w+1
!            ww(counter) = tw(i+kw-1) + (real(j)/(n_sub_w+1)*(tw(i+kw)-tw(i+kw-1)))
!            counter = counter + 1
!         end do
!      end do

!      allocate(volume_vals(ndim,nww,nvv,nuu))

!      do i=1,nuu
!         do j=1,nvv
!            do k=1,nww
!               call eval_volume(uu(i),vv(j),ww(k),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,&
!                    1,1,1,volume_vals(:,k,j,i))
!            end do
!         end do
!      end do
!      ! Now do comparison
!      u0(:) = u(:)   
!      v0(:) = v(:)
!      w0(:) = w(:)
!      do ipt=1,N
!         if (u(ipt) < 0 .or. u(ipt) > 1 .or. v(ipt) < 0 .or. v(ipt) > 1 .or. w(ipt) < 0 .or. w(ipt) > 1) then
!            val0 = volume_vals(:,1,1,1)
!            D0 = norm(val0-x0(ipt,:),ndim)
!            u0(ipt) = uu(1)
!            v0(ipt) = vv(1)
!            w0(ipt) = ww(1)
!            do i = 1,nuu
!               do j = 1,nvv
!                  do k = 1,nww

!                     D = norm(volume_vals(:,k,j,i)-X0(:,ipt),ndim)
!                     if (D<D0) then
!                        u0(ipt) = uu(i)
!                        v0(ipt) = vv(j)
!                        w0(ipt) = ww(k)
!                        D0 = D
!                     end if
!                  end do
!               end do
!            end do
!         else
!            u0(ipt) = u(ipt)
!            v0(ipt) = v(ipt)
!            w0(ipt) = w(ipt)
!         end if
!      end do
!   else
!      u0(:) = u(:)   
!      v0(:) = v(:)
!      w0(:) = w(:)
!   end if

!   u = u0(:)
!   v = v0(:)
!   w = w0(:)

!   ! Now we have u0,v0,w0 so we can do the newton search
!   do ipt=1,N

!      call eval_volume(u(ipt),v(ipt),w(ipt),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,1,1,1,val)
!      call eval_volume_deriv(u(ipt),v(ipt),w(ipt),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,deriv)
!      call eval_volume_deriv2(u(ipt),v(ipt),w(ipt),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,deriv2)

!      Diff(:,ipt) = val-x0(:,ipt)

!      iteration_loop: do i=1,niter
!         ! Check the convergence criteria
!         if (norm(Diff(:,ipt),ndim) <= eps1) then
           
!            exit iteration_loop
!         end if

!         c1 = ndp(deriv(:,1),diff(:,ipt),ndim)/(norm(deriv(:,1),ndim)*norm(Diff(:,ipt),ndim))
!         c2 = ndp(deriv(:,2),diff(:,ipt),ndim)/(norm(deriv(:,2),ndim)*norm(Diff(:,ipt),ndim))
!         c3 = ndp(deriv(:,3),diff(:,ipt),ndim)/(norm(deriv(:,3),ndim)*norm(Diff(:,ipt),ndim))
 
!         if (c1 <= eps2 .and. c2 <= eps2 .and. c3 <=eps2) then
!            exit iteration_loop
!         end if
!         u0(ipt) = u(ipt)
!         v0(ipt) = v(ipt)
!         w0(ipt) = w(ipt)
!         call eval_volume(u(ipt),v(ipt),w(ipt),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,1,1,1,val)
!         call eval_volume_deriv(u(ipt),v(ipt),w(ipt),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,deriv)
!         call eval_volume_deriv2(u(ipt),v(ipt),w(ipt),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,deriv2)

!         Diff(:,ipt) = val-x0(:,ipt)

!         A(1,1) = norm(deriv(:,1),ndim)**2 + dotproduct(Diff(:,ipt),deriv2(:,1,1),3)
!         A(1,2) = dotproduct(deriv(:,1),deriv(:,2),3) + dotproduct(Diff(:,ipt),deriv2(:,1,2),3)
!         A(1,3) = dotproduct(deriv(:,1),deriv(:,3),3) + dotproduct(Diff(:,ipt),deriv2(:,1,3),3)
!         A(2,2) = norm(deriv(:,2),ndim)**2 + dotproduct(Diff(:,ipt),deriv2(:,2,2),3)
!         A(2,3) = dotproduct(deriv(:,2),deriv(:,3),3) + dotproduct(Diff(:,ipt),deriv2(:,2,3),3)
!         A(3,3) = norm(deriv(:,3),ndim)**2 + dotproduct(Diff(:,ipt),deriv2(:,3,3),3)
!         A(2,1) = A(1,2)
!         A(3,1) = A(1,3)
!         A(3,2) = A(2,3)

!         ki(1) = -dotproduct(Diff(:,ipt),deriv(:,1),3)
!         ki(2) = -dotproduct(Diff(:,ipt),deriv(:,2),3)
!         ki(3) = -dotproduct(Diff(:,ipt),deriv(:,3),3)
      
!         call solve_3by3(A,ki,delta)

!         ! -- U Bounds Checking --
!         if (u0(ipt)+delta(1)  < tu(1)) then
!            delta(1) = tu(1)-u0(ipt)
!         end if

!         if (u0(ipt)+delta(1) > tu(nctlu+ku)) then
!            delta(1) = tu(nctlu+ku)-u0(ipt)
!         end if

!         ! -- V Bounds Checking --
!         if (v0(ipt)+delta(2) < tv(1)) then
!            delta(2) = tv(1)-v0(ipt)
!         end if

!         if (v0(ipt)+delta(2) > tv(nctlv+kv)) then
!            delta(2) = tv(nctlv+kv)-v0(ipt)
!         end if
        
!         ! -- W Bounds Checking --
!         if (w0(ipt)+delta(3) < tw(1)) then
!            delta(3) = tw(1)-w0(ipt)
!         end if

!         if (w0(ipt)+delta(3) > tw(nctlw+kw)) then
!            delta(3) = tw(nctlw+kw)-w0(ipt)
!         end if

!         inner_loop: do j=1,max_inner_iter
!            u(ipt) = u0(ipt) + delta(1)
!            v(ipt) = v0(ipt) + delta(2)
!            w(ipt) = w0(ipt) + delta(3)
           
!            call eval_volume(u(ipt),v(ipt),w(ipt),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,1,1,1,val)

!            D2 = val-x0(:,ipt)

!            if ( (norm(D2,ndim) -norm(diff(:,ipt),ndim)) .ge. eps1) then
!               delta = delta * 0.5
!            else
!               exit inner_loop
!            end if
!         end do inner_loop

!         ! No Change convergence Test

!         if (norm( (u(ipt)-u0(ipt))*deriv(:,1) + &
!              (v(ipt)-v0(ipt))*deriv(:,2) + &
!              (w(ipt)-w0(ipt))*deriv(:,3),ndim) <= eps1) then
!            exit iteration_loop
!         end if
!      end do iteration_loop
!   end do
!   if (brute_force) then
!      deallocate(volume_vals)
!      deallocate(uu,vv,ww)
!   end if
 
! end subroutine point_volume_old

subroutine curve_curve(t1,k1,coef1,t2,k2,coef2,n1,n2,ndim,Niter,eps1,eps2,s,t,Diff)

  !*** DESCRIPTION
  !
  !     Written by Gaetan Kenway
  ! 
  !     The function takes the spline defination of a bivariate spline
  !     surface as defined by coef,kx,ky,nx,ny,tx,ty and a point x0 and
  !     determines the u,v positions that minimizes the distance between
  !     the point and the surface. 

  !     Description of Arguments:
  !     Input:
  !     t1      - Knot Vector for Curve 1
  !     k1      - Order for Curve 1
  !     coef1   - Coefficients for Curve 1
  !     t2      - Knot Vector for Curve 2
  !     k2      - Order for Curve 2
  !     coef2   - Coefficients for Curve 2
  !     n1      - Number of coefficients on curve 1
  !     n2      - Number of coefficients on curve 2
  !     Niter   - Integer: Maximum number of iterations
  !     tol     - Real: Tolerance for newton iteration
  !     s       - Real: Initial guess for parameter on Curve 1
  !     t       - Real: Initial guess for parameter on Curve 2

  !     Output:
  !     s       - Real: parameter on Curve 1
  !     t       - Real: parameter on Curve 2
  !     D       - Real: Distance between curve 
  !
  implicit none

  ! Input/Output
  integer         , intent(in)     :: n1,n2,k1,k2,ndim
  double precision, intent(in)     :: t1(n1+k1),t2(n2+k2)
  double precision, intent(in)     :: coef1(ndim,n1),coef2(ndim,n2)
  double precision, intent(inout)  :: s,t 
  integer         , intent(in)     :: Niter
  double precision, intent(in)     :: eps1,eps2

  double precision, intent(out)    :: Diff(ndim)

  ! Working
  double precision                 :: val0(ndim),val1(ndim),val2(ndim)
  double precision                 :: deriv_c1(ndim),deriv2_c1(ndim)
  double precision                 :: deriv_c2(ndim),deriv2_c2(ndim)
  integer                          :: i,j,ii,jj,max_inner_iter,counter1,istart
  double precision                 :: D,D0,s0,t0,delta(2),D2(3)
  double precision                 :: ki(2),A(2,2),val0_1(3),val0_2(3)
  integer                          :: n,nss,ntt ! Huristic Value
  double precision,allocatable     :: curve1_vals(:,:),curve2_vals(:,:)
  double precision,allocatable     :: ss(:),tt(:)
  logical                          :: brute_force
  ! Functions 
  double precision                 :: norm
  
  max_inner_iter = 20
  ! Is 3 Good Here?? Huristic
  if (k1 == 2 .or. k2 == 2) then
     n = 100
  else
     n = 100
  end if
  brute_force = .False.
  if (s < 0 .or. s > 1 .or. t < 0 .or. t > 1) then
     ! Do a brute force approach to get good starting point
     brute_force = .True.
     nss = n1-k1+2 + (n1-k1+1)*n
     ntt = n2-k2+2 + (n2-k2+1)*n
     allocate (ss(nss),tt(ntt))

     counter1 = 1
     do i=1,n1-k1+1
        if (i==1) then
           istart = 0
        else
           istart = 1
        end if
        do j=istart,n+1
           ss(counter1) = t1(i+k1-1) + (real(j)/(n+1)*(t1(i+k1)-t1(i+k1-1)))
           counter1 = counter1 + 1
        end do
     end do

     counter1 = 1
     do i=1,n2-k2+1
        if (i==1) then
           istart = 0
        else
           istart = 1
        end if
        do j=istart,n+1
           tt(counter1) = t2(i+k2-1) + (real(j)/(n+1)*(t2(i+k2)-t2(i+k2-1)))
           counter1 = counter1 + 1
        end do
     end do

     allocate(curve1_vals(ndim,nss),curve2_vals(ndim,ntt))
     do i=1,nss
        call eval_curve(ss(i),t1,k1,coef1,n1,ndim,1,curve1_vals(:,i))
     end do

     do i=1,ntt
        call eval_curve(tt(i),t2,k2,coef2,n2,ndim,1,curve2_vals(:,i))
     end do
     
     val0_1 = curve1_vals(:,1)
     val0_2 = curve2_vals(:,1)
     D0 = norm(val0_1-val0_2,ndim)
     s0 = ss(1)
     t0 = tt(1)
     do i = 1,nss
        do j = 1,ntt
           D = norm(curve1_vals(:,i)-curve2_vals(:,j),ndim)
           if (D<D0) then
              s0 = ss(i)
              t0 = tt(j)
              D0 = D
           end if
        end do
     end do
  else
     s0 = s
     t0 = t
     call eval_curve(s,t1,k1,coef1,n1,ndim,1,val1)
     call eval_curve(t,t2,k2,coef2,n2,ndim,1,val2)
     D0 = norm(val1-val2,ndim)
  end if

  ! Now we have s0 and t0 so we can do the newton iteration
  call eval_curve(s0,t1,k1,coef1,n1,ndim,1,val1)
  call eval_curve_deriv(s0,t1,k2,coef1,n1,ndim,deriv_c1)
  call eval_curve_deriv2(s0,t1,k2,coef1,n1,ndim,deriv2_c1)

  call eval_curve(t0,t2,k2,coef2,n2,ndim,1,val2)
  call eval_curve_deriv(t0,t2,k2,coef2,n2,ndim,deriv_c2)
  call eval_curve_deriv2(t0,t2,k2,coef2,n2,ndim,deriv2_c2)
  Diff = val1-val2
  s = s0
  t = t0

  iteration_loop: do i=1,niter
     ! Check the converge criteria -> coincidence
     if (norm(Diff,ndim) <= eps1) then
        exit iteration_loop
     end if

     ! Cosine convergence check

     if ( norm(dot_product(deriv_c1,Diff),ndim)/(norm(deriv_c1,ndim)*norm(Diff,ndim)) < eps2 .and. &
          norm(dot_product(deriv_c2,Diff),ndim)/(norm(deriv_c2,ndim)*norm(Diff,ndim)) < eps2) then
        exit iteration_loop
     end if
     s0 = s
     t0 = t

     call eval_curve(s,t1,k1,coef1,n1,ndim,1,val1)
     call eval_curve_deriv(s,t1,k2,coef1,n1,ndim,deriv_c1)
     call eval_curve_deriv2(s,t1,k2,coef1,n1,ndim,deriv2_c1)

     call eval_curve(t,t2,k2,coef2,n2,ndim,1,val2)
     call eval_curve_deriv(t,t2,k2,coef2,n2,ndim,deriv_c2)
     call eval_curve_deriv2(t,t2,k2,coef2,n2,ndim,deriv2_c2)

     Diff = val1-val2

     A(1,1) = norm(deriv_c1,ndim)**2 + dot_product(Diff,deriv2_c1)
     A(2,2) = -dot_product(deriv_c1,deriv_c2)
     A(2,1) = dot_product(deriv_c1,deriv_c2)
     A(2,2) = -norm(deriv_c2,ndim)**2 + dot_product(Diff,deriv2_c2)

     ki(1) = -dot_product(Diff,deriv_c1)
     ki(2) = -dot_product(Diff,deriv_c2)

     call solve_2by2(A,ki,delta)

     ! Bounds checking
     if (s+delta(1) < t1(1)) then
        delta(1) = t1(1)-s
     end if
     if (s+delta(1) > t1(n1+k1)) then
        delta(1) = t1(n1+k1) - s
     end if

     if (t+delta(2) < t2(1)) then
        delta(2) = t2(1)-t
     end if
     if (t+delta(2) > t2(n2+k2)) then
        delta(2) = t2(n2+k2) - t
     end if

     inner_loop: do j=1,max_inner_iter
        s = s0 + delta(1)
        t = t0 + delta(2)
        call eval_curve(s,t1,k1,coef1,n1,ndim,1,val1)
        call eval_curve(t,t2,k2,coef2,n2,ndim,1,val2)
        D2 = val1-val2
        if ((norm(D2,ndim)-norm(Diff,ndim)) .ge. eps1) then
           delta = 0.5*delta 
        else
           exit inner_loop
        end if
     end do inner_loop

     ! No change convergence Test
     if (norm( (s-s0)*deriv_c1 + (t-t0)*deriv_c2,ndim) <= eps1) then
        exit iteration_loop
     end if

  end do iteration_loop

  if (brute_force) then
     deallocate(curve1_vals,curve2_vals,ss,tt)
  end if

end subroutine curve_curve

subroutine curve_surface(tc,kc,coefc,tu,tv,ku,kv,coefs,nctlc,nctlu,nctlv,ndim,niter,eps1,eps2,u,v,s,Diff)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: curve_surface attempts to find the minimum distance between 
  !     a curve and a surface. This works for intersections as well
  !
  !     Description of Arguments
  !     Input
  !     ct      - Real, curve knot vector. Length nctlc+kc
  !     kc      - Integer, curve order
  !     coefc   - Real, aray size(ndim,nctlc) coefficients for curve
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     coefs   - Real,Array of B-spline coefficients  Size (ndim,nctlv,nctlu)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     ndim    - Integer, Spatial Dimension
  !     Niter   - Integer, Maximum number of Netwton iterations
  !     eps1    - Real - Eculdian Distance Convergence Measure
  !     eps2    - Real - Cosine Convergence Measure
  !
  !     Ouput 
  !     u       - Real, u parameter where S(u,v) is closest to C(s)
  !     v       - Real, v parameter where S(u,v) is closest to C(s)
  !     s       - Real, s parameter where S(u,v) is closest to C(s)
  !     diff    - Real Array size(ndim) - Distance between C(s)and S(u,v)

  implicit none
  ! Input
  integer         , intent(in)          :: kc,ku,kv,nctlu,nctlv,nctlc,ndim,niter
  double precision, intent(in)          :: tc(nctlc+kc),tu(nctlu+ku),tv(nctlv+kv)
  double precision, intent(in)          :: coefc(ndim,nctlc),coefs(ndim,nctlv,nctlu)
  double precision, intent(in)          :: eps1,eps2
  ! Output
  double precision, intent(inout)       :: u,v,s
  double precision, intent(out)         :: diff(ndim)

  ! Working
  double precision                      :: val_s(ndim),deriv_s(ndim,2),deriv2_s(ndim,2,2)
  double precision                      :: val_c(ndim),deriv_c(ndim),deriv2_C(ndim)
  double precision                      :: val0_s(ndim),val0_c(ndim)
  integer                               :: i,j,ii,jj,l,ll,counter1,counter2
  double precision                      :: D,D0,u0,v0,s0,delta(3),D2(3)
  double precision                      :: A(3,3),ki(3)
  integer                               :: n_sub_u,n_sub_v,nc ! Huristic Value
  integer                               :: istart,nss,nuu,nvv,max_inner_iter
  logical                               :: brute_force
  ! Allocatable 
  double precision,allocatable          :: curve_vals(:,:),uu(:),vv(:),ss(:)
  double precision,allocatable          :: surface_vals(:,:,:)

  ! Functions     
  double precision                      :: norm

 max_inner_iter = 20

 if (ku == 2) then
     n_sub_u = 3
  else
     n_sub_u = 3
  end if

  if (kv == 2) then
     n_sub_v = 3
  else
     n_sub_v = 3
  end if

  if (kc == 2) then
     nc = 3
  else
     nc = 10
  end if

  ! First we will evaluate the surface at n points inside each knot span in each direction
  !if we are given a bad guess do the brute force
  brute_force = .False.
  if (u < 0 .or. u > 1 .or. v < 0 .or. v > 1 .or. s < 0 .or. s > 1) then
     brute_force = .True.

     ! Generate the ss,uu and vv values
     nss = nctlc-kc+2 + (nctlc-kc+1)*nc
     nuu = nctlu-ku+2 + (nctlu-ku+1)*n_sub_u
     nvv = nctlv-kv+2 + (nctlv-kv+1)*n_sub_v
     allocate (ss(nss),uu(nuu),vv(nvv))

     counter1 = 1
     do i=1,nctlc-kc+1
        if (i==1) then
           istart = 0
        else
           istart = 1
        end if
        do j=istart,nc+1
           ss(counter1) = tc(i+kc-1) + (real(j)/(nc+1)*(tc(i+kc)-tc(i+kc-1)))
           counter1 = counter1 + 1
        end do
     end do

     counter1 = 1
     do i=1,nctlu-ku+1
        if (i==1) then
           istart = 0
        else
           istart = 1
        end if
        do j=istart,n_sub_u+1
           uu(counter1) = tu(i+ku-1) + (real(j)/(n_sub_u+1)*(tu(i+ku)-tu(i+ku-1)))
           counter1 = counter1 + 1
        end do
     end do

     counter1 = 1
     do i=1,nctlv-kv+1
        if (i==1) then
           istart = 0
        else
           istart = 1
        end if
        do j=istart,n_sub_v+1
           vv(counter1) = tv(i+kv-1) + (real(j)/(n_sub_v+1)*(tv(i+kv)-tv(i+kv-1)))
           counter1 = counter1 + 1
        end do
     end do

     allocate(curve_vals(ndim,nss),surface_vals(ndim,nvv,nuu))

     do i=1,nss
        call eval_curve(ss(i),tc,kc,coefc,nctlc,ndim,1,curve_vals(:,i))
     end  do

     do i=1,nuu
        do j=1,nvv
           call eval_surface(uu(i),vv(j),tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,1,1,&
                surface_vals(:,j,i))
        end do
     end do

     ! Now do comparison
     val0_s = surface_vals(:,1,1)
     val0_c = curve_vals(:,1)
     D0 = norm(val0_s-val0_c,ndim)
     u0 = uu(1)
     v0 = vv(1)
     s0 = ss(1)
     do l = 1,nss
        do i = 1,nuu
           do j = 1,nvv
              D = norm(surface_vals(:,j,i)-curve_vals(:,l),ndim)
              if (D<D0) then
                 u0 = uu(i)
                 v0 = vv(j)
                 s0 = ss(l)
                 D0 = D
              end if
           end do
        end do
     end do
  else
     u0 = u
     v0 = v
     s0 = s
  end if

  ! Now we have u0,v0,s0 so we can do the newton iteration

  call eval_surface(u0,v0,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,1,1,val_s)
  call eval_surface_deriv(u0,v0,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,deriv_s)
  call eval_surface_deriv2(u0,v0,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,deriv2_s)

  call eval_curve(s0,tc,kc,coefc,nctlc,ndim,1,val_c)
  call eval_curve_deriv(s0,tc,kc,coefc,nctlc,ndim,deriv_c)
  call eval_curve_deriv2(s0,tc,kc,coefc,nctlc,ndim,deriv2_c)
  
  Diff = val_s-val_c

  u = u0
  v = v0
  s = s0

  iteration_loop: do i=1,niter

     ! Check the convergence criteria
     if (norm(Diff,ndim) <= eps1) then
        exit iteration_loop
     end if

     if ( norm(dot_product(deriv_s(:,1),Diff),ndim)/(norm(deriv_s(:,1),ndim)*norm(Diff,ndim)) <= eps2 .and. &
          norm(dot_product(deriv_s(:,2),Diff),ndim)/(norm(deriv_s(:,2),ndim)*norm(Diff,ndim)) <= eps2 .and. &
          norm(dot_product(deriv_c     ,Diff),ndim)/(norm(deriv_c     ,ndim)*norm(Diff,ndim)) <= eps2) then
        exit iteration_loop
     end if
     u0 = u
     v0 = v
     s0 = s

     call eval_surface(u0,v0,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,1,1,val_s)
     call eval_surface_deriv(u0,v0,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,deriv_s)
     call eval_surface_deriv2(u0,v0,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,deriv2_s)

     call eval_curve(s0,tc,kc,coefc,nctlc,ndim,1,val_c)
     call eval_curve_deriv(s0,tc,kc,coefc,nctlc,ndim,deriv_c)
     call eval_curve_deriv2(s0,tc,kc,coefc,nctlc,ndim,deriv2_c)

     Diff = val_s-val_c

     A(1,1) = norm(deriv_s(:,1),ndim)**2 + dot_product(Diff,deriv2_s(:,1,1))
     A(1,2) = dot_product(deriv_s(:,1),deriv_s(:,2)) + dot_product(Diff,deriv2_s(:,1,2))
     A(1,3) = -dot_product(deriv_c,deriv_s(:,1))

     A(2,1) = A(1,2)
     A(2,2) = norm(deriv_s(:,2),ndim)**2 + dot_product(Diff,deriv2_s(:,2,2))
     A(2,3) = -dot_product(deriv_c,deriv_s(:,2))

     A(3,1) = dot_product(deriv_s(:,1),deriv_c) + dot_product(Diff,deriv2_c)
     A(3,2) = dot_product(deriv_s(:,2),deriv_c) + dot_product(Diff,deriv2_c)
     A(3,3) = -norm(deriv_c,ndim)**2 + dot_product(Diff,deriv2_c)

     ki(1) = -dot_product(Diff,deriv_s(:,1))
     ki(2) = -dot_product(Diff,deriv_s(:,2))
     ki(3) = -dot_product(Diff,deriv_c)
     call solve_3by3(A,ki,delta)

     ! Error Checking

     ! Bounds Checking
     if (u0+delta(1) < tu(1)) then
        delta(1) = tu(1)-u0
     end if

     if (u0+delta(1) > tu(nctlu+ku)) then
        delta(1) = tu(nctlu+ku)-u0
     end if

     if (v0+delta(2) < tv(1)) then
        delta(2) = tv(1)-v0
     end if

     if (v0+delta(2) > tv(nctlv+kv)) then
        delta(2) = tv(nctlv+kv)-v0
     end if

     if (s0+delta(3) < tc(1)) then
        delta(3) = tc(1)-s0
     end if

     if (s0+delta(3) > tc(nctlc+kc)) then
        delta(3) = tc(nctlc+kc) - s0
     end if

     inner_loop: do j=1,max_inner_iter
        u = u0 + delta(1)
        v = v0 + delta(2)
        s = s0 + delta(3)
        call eval_surface(u,v,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,1,1,val_s)
        call eval_curve(s,tc,kc,coefc,nctlc,ndim,1,val_c)

        D2 = val_s-val_c
        if ((norm(D2,ndim)-norm(Diff,ndim)) .ge. eps1) then
           delta = 0.5*delta 
        else
           exit inner_loop
        end if
     end do inner_loop
        
     ! No Change convergence Test

     if (norm((u-u0)*deriv_s(:,1)+(v-v0)*deriv_s(:,2)+(s-s0)*deriv_c,ndim) <= eps1) then
        exit iteration_loop
     end if

  end do iteration_loop

  if (brute_force) then
     deallocate(curve_vals,surface_vals,ss,uu,vv)
  end if

end subroutine curve_surface

subroutine solve_2by2(A,b,x)
  ! Solve a 2 x 2 system  -- With NO checking
  double precision, intent(in) :: A(2,2),b(2)
  double precision, intent(out) :: x(2)
  double precision         :: idet

  idet = 1/(A(1,1)*A(2,2)-A(1,2)*A(2,1))
  x(1) = idet*(A(2,2)*b(1) - A(1,2)*b(2))
  x(2) = idet*(-A(2,1)*b(1) + A(1,1)*b(2))

end subroutine solve_2by2

subroutine solve_3by3(A,b,x)
  ! Solve a 3 x 3 system  -- With NO checking
  double precision,intent(in) :: A(3,3),b(3)
  double precision,intent(out) :: x(3)
  double precision         :: idet

  idet = 1/(A(1,1)*(A(3,3)*A(2,2)-A(3,2)*A(2,3))-A(2,1)*(A(3,3)*A(1,2)-A(3,2)*A(1,3))+A(3,1)*(A(2,3)*A(1,2)-A(2,2)*A(1,3)))
  x(1) = idet*( b(1)*(A(3,3)*A(2,2)-A(3,2)*A(2,3)) - b(2)*(A(3,3)*A(1,2)-A(3,2)*A(1,3)) + b(3)*(A(2,3)*A(1,2)-A(2,2)*A(1,3)))
  x(2) = idet*(-b(1)*(A(3,3)*A(2,1)-A(3,1)*A(2,3)) + b(2)*(A(3,3)*A(1,1)-A(3,1)*A(1,3)) - b(3)*(A(2,3)*A(1,1)-A(2,1)*A(1,3)))
  x(3) = idet*( b(1)*(A(3,2)*A(2,1)-A(3,1)*A(2,2)) - b(2)*(A(3,2)*A(1,1)-A(3,1)*A(1,2)) + b(3)*(A(2,2)*A(1,1)-A(2,1)*A(1,2)))


  ! | a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
  ! | a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
  ! | a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |

  ! DET  =  a11(a33a22-a32a23)-a21(a33a12-a32a13)+a31(a23a12-a22a13)

end subroutine solve_3by3

subroutine line_plane(ia,vc,p0,v1,v2,n,sol,n_sol)

  ! Check a line against multiple planes

  implicit none 
  ! Input
  integer, intent(in) :: n
  double precision , intent(in) :: ia(3),vc(3),p0(3,n),v1(3,n),v2(3,n)

  ! Output
  integer, intent(out) :: n_sol
  double precision, intent(out) :: sol(6,n)

  ! Worling 
  integer :: i,ind(n)
  double precision :: A(3,3),rhs(3),x(3),norm
  
  A(:,1) = -vc
  n_sol = 0
  sol(:,:) = 0.0
  do i=1,n
     A(:,2) = v1(:,i)
     A(:,3) = v2(:,i)
     rhs = ia-p0(:,i)
     
     call solve_3by3(A,rhs,x)
     
     if (x(1) .ge. 0.00 .and. x(1) .le. 1.00 .and. &
         x(2) .ge. 0.00 .and. x(2) .le. 1.00 .and. &
         x(3) .ge. 0.00 .and. x(3) .le. 1.00 .and. &
         x(2)+x(3) .le. 1.00) then

        n_sol = n_sol + 1
        sol(1:3,n_sol) = x  ! t,u,v parametric locations
        sol(4:6,n_sol) = ia + x(1)*vc ! Actual point value
        ind(n_sol) = i

     end if
  end do
end subroutine line_plane


subroutine plane_line(ia,vc,p0,v1,v2,n,sol,n_sol)

  ! Check a plane against multiple lines

  implicit none 
  ! Input
  integer, intent(in) :: n
  double precision , intent(in) :: ia(3,n),vc(3,n),p0(3),v1(3),v2(3)

  ! Output
  integer, intent(out) :: n_sol
  double precision, intent(out) :: sol(6,n)

  ! Worling 
  integer :: i,ind(n)
  double precision :: A(3,3),rhs(3),x(3),norm

  n_sol = 0
  sol(:,:) = 0.0
  A(:,2) = v1(:)
  A(:,3) = v2(:)
  
  do i=1,n
  
     A(:,1) = -vc(:,i)
     rhs = ia(:,i)-p0(:)
     
     call solve_3by3(A,rhs,x)
     
     if (x(1) .ge. 0.00 .and. x(1) .le. 1.00 .and. &
         x(2) .ge. 0.00 .and. x(2) .le. 1.00 .and. &
         x(3) .ge. 0.00 .and. x(3) .le. 1.00 .and. &
         x(2)+x(3) .le. 1.00) then

        n_sol = n_sol + 1
        sol(1:3,n_sol) = x  ! t,u,v parametric locations
        sol(4:6,n_sol) = ia(:,i) + x(1)*vc(:,i) ! Actual point value
        ind(n_sol) = i

     end if
  end do
end subroutine plane_line

! subroutine find_volume(p0,v1,v2,x0,lines,n,nVolids)
!   ! Determine what volume defined by 12 triangular faces in p0,v1,v2
!   ! the points x0 coorespond to

!   implicit none 
!   ! Input
!   integer, intent(in) :: n,nVol
!   double precision , intent(in) :: p0(3,12*nVol),v1(3,12*nVol),v2(3,12*nVol)
!   double precision , intent(in) :: x0(3,n)
!   double precision , intent(in) :: lines(3,6)

!   ! Output
!   integer, intent(out) :: ids(n)

!   ! Working 
!   integer :: i,ind(n)
!   double precision :: A(3,3),rhs(3),x(3),norm
!   double precision :: sol(12)
!   double precision :: D(n)

!   ! Loop over each point:
!   do i=1,n
     
!      vol_loop: do iVol=1,nVol

!         ! Loop over each face
!         do j=1,12

!            A(:,2) = v1(:,(iVol-1)*12+j)
!            A(:,3) = v2(:,(iVol-1)*12+j)

!            ! Loop over each line:
!            do k=1,6
                       
!               A(:,1) = -lines(:,k)

!               rhs = x0(:,i)-p0(:)
        
!               call solve_3by3(A,rhs,x)
        
!               if (x(1) .ge. 0.00 .and. x(1) .le. 1.00 .and. &
!                    x(2) .ge. 0.00 .and. x(2) .le. 1.00 .and. &
!                    x(3) .ge. 0.00 .and. x(3) .le. 1.00 .and. &
!                    x(2)+x(3) .le. 1.00) then
!                  sol(j) = 1
              
!                  ! Distance to patch is x(1) = t
!               end if
!            end do
!         if (sum(sol) >= 6) then
!            ! We've definately found the right volume...set it:
!            ids(i) = iVol
!            exit vol_loop
!         end if

!      end do vol_loop



              

! !        sol(1:3,n_sol) = x  ! t,u,v parametric locations
! !        sol(4:6,n_sol) = ia(:,i) + x(1)*vc(:,i) ! Actual point value
! !        ind(n_sol) = i

!      end if
!   end do



! end subroutine find_volume

  
subroutine point_plane(pt,p0,v1,v2,n,sol,n_sol,best_sol)

  implicit none 
  ! Input
  integer, intent(in) :: n
  double precision , intent(in) :: pt(3),p0(3,n),v1(3,n),v2(3,n)

  ! Output
  integer, intent(out) :: n_sol,best_sol
  double precision, intent(out) :: sol(6,n)

  ! Worling 
  integer :: i,ind(n)
  double precision :: A(2,2),rhs(2),x(2),norm,r(3),D,D0
  

  n_sol = 0
  sol(:,:) = 0.0
  do i=1,n
     A(1,1) = v1(1,i)**2 + v1(2,i)**2 + v1(3,i)**2
     A(1,2) = v1(1,i)*v2(1,i) + v1(2,i)*v2(2,i) + v1(3,i)+v2(3,i)
     A(2,1) = A(1,2)
     A(2,2) = v2(1,i)**2 + v2(2,i)**2 + v2(3,i)**2
     r = p0(:,i)-pt
     rhs(1) = r(1)*v1(1,i) + r(2)*v1(2,i) + r(3)*v1(3,i)
     rhs(2) = r(1)*v2(1,i) + r(2)*v2(2,i) + r(3)*v2(3,i)
     
     call solve_2by2(A,rhs,x)
     
     if (x(1) .ge. 0.00 .and. x(1) .le. 1.00 .and. &
         x(2) .ge. 0.00 .and. x(2) .le. 1.00 .and. &
         x(1)+x(2) .le. 1.00) then

        n_sol = n_sol + 1
        sol(2:3,n_sol) = x  ! t,u,v parametric locations
        sol(4:6,n_sol) = 0.0 ! Actual point value
        ind(n_sol) = i
     end if
  end do

  ! Now post-process to get the closest one
  best_sol = 1
  D0 = norm(p0(:,ind(1)) + sol(2,ind(1))*v1(:,ind(1)) + sol(3,ind(1))*v2(:,ind(1)))

  do i=1,n_sol
     D = norm(p0(:,ind(i)) + sol(2,ind(i))*v1(:,ind(i)) + sol(3,ind(i))*v2(:,ind(i)),3)
     if (D<D0) then
        D0 = D
        best_sol = i
     end if
  end do
end subroutine point_plane
     
  
function dotproduct(x1,x2,n)
  implicit none
  double precision,intent(in) :: x1(n),x2(n)
  integer,intent(in) :: n
  integer :: i

  double precision ::  dotproduct 
  dotproduct = 0.0
  do i=1,n
     dotproduct = dotproduct + x1(i)*x2(i)
  end do

end function dotproduct
  
function ndp(x1,x2,n)
  ! Compute norm of the dot_product
  implicit none
  double precision,intent(in) :: x1(n),x2(n)
  integer,intent(in) :: n
  integer :: i

  double precision ::  ndp
  ndp = 0.0
  do i=1,n
     ndp = ndp + (x1(i)*x2(i))**2
  end do
  ndp = sqrt(ndp)
end function ndp
  
