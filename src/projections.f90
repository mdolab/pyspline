subroutine point_curve(x0,t,k,coef,nctl,ndim,N,Niter,eps1,eps2,s,Diff)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: point_curve attempts to solve the point inversion problem
  !
  !     Description of Arguments
  !     Input
  !     x0      - Real, array size(N,ndim) N points we are trying to invert
  !     t       - Real, Knot vector. Length nctl+k
  !     k       - Integer,order of B-spline
  !     coef    - Real,Array of B-spline coefficients and weights. Size (nctl,ndim)
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
  !     diff    - Real,array, size(N,ndim)- Distance between x0 and curve(s)

  implicit none
  ! Input
  double precision, intent(in)          :: x0(N,ndim)
  integer         , intent(in)          :: k,nctl,ndim,N,niter
  double precision, intent(in)          :: t(nctl+k)
  double precision, intent(in)          :: coef(nctl,ndim)
  double precision, intent(in)          :: eps1,eps2

  ! Output
  double precision, intent(out)         :: s(N),Diff(N,ndim)

  ! Working
  double precision                      :: val(ndim),deriv(ndim),deriv2(ndim)
  double precision                      :: ss
  double precision                      :: val0(ndim),s0(N)
  integer                               :: i,j,ipt,counter
  double precision                      :: D,D0
  integer                               :: n_sub ! Huristic Value
  logical                               :: brute_force
  ! Alloctable
  double precision,allocatable          :: curve_vals(:,:)

  ! Functions
  double precision                      :: norm

  ! Initialization
  n_sub=10  ! Is 3  Good Here?? Huristic!
  brute_force = .False.
  ! Determine if any of the guesses are out of range, if so, do the brute force search

  point_loop: do ipt=1,N
     if (s(ipt) < 0 .or. s(ipt) > 1) then
        brute_force = .True.
        exit point_loop
     end if
  end do point_loop

  if (brute_force .eqv. .True.) then
     ! Dynamically allot the required space for the brute-force search search
     allocate(curve_vals((nctl-k+2) + (nctl-k+1)*n_sub,ndim))
     counter = 1
     do i=1,nctl-k+2 ! Number of knot spans for k-repeated knots at end
        call eval_curve(t(i+k-1),t,k,coef,nctl,ndim,curve_vals(counter,:))
        counter = counter + 1
     end do
     do i=1,nctl-k+1
        do j =1,n_sub
           ss = t(i+k-1) + (real(j)/(n_sub+1)*(t(i+k)-t(i+k-1)))
           counter = counter + 1
        end do
     end do
     s0(:) = s(:)
     ! Now do the actual searching
     do ipt=1,N
        if (s(ipt) < 0 .or. s(ipt) > 1) then ! Still only do it if we have a bad guess
           call eval_curve(t(1),t,k,coef,nctl,ndim,val0)
           D0 = norm(val0-x0(ipt,:),ndim)
           counter = 1
           do i=1,nctl-k+2
              D  = norm(x0(ipt,:)-curve_vals(counter,:),ndim)
              counter = counter + 1
              if (D<D0) then
                 s0(ipt) = t(i+k-1)
                 D0 = D
              end if
           end do

           do i=1,nctl-k+1
              do j=1,n_sub
                 D  = norm(x0(ipt,:)-curve_vals(counter,:),ndim)
                 counter = counter + 1
                 if (D<D0) then
                    s0(ipt) = t(i+k-1) + (real(j)/(n_sub+1)*(t(i+k)-t(i+k-1)))
                    D0 = D
                 end if
              end do
           end do
        else
           s0(i) = s(i)
        end if
     end do
  else
     s0(:) = s(:)   
  end if

  do ipt=1,N
     ! Now we have s0 so we should be able to do the netwon seach
     call eval_curve(s0(ipt),t,k,coef,nctl,ndim,val)
     call eval_curve_deriv(s0(ipt),t,k,coef,nctl,ndim,deriv)
     call eval_curve_deriv2(s0(ipt),t,k,coef,nctl,ndim,deriv2)
     Diff(ipt,:) = val-x0(ipt,:)
     s(ipt) = s0(ipt) - dot_product(deriv,Diff(ipt,:))/(dot_product(deriv2,Diff(ipt,:)) + norm(deriv,ndim)**2)
     iteration_loop: do i=1,Niter
        ! Check the Convergence Criteria
        if (norm(Diff(ipt,:),ndim) <= eps1) then
           exit iteration_loop
        end if

        if (  norm(dot_product(deriv,Diff(ipt,:)),ndim)/ (norm(deriv,ndim)*norm(Diff(ipt,:),ndim)) <= eps2) then
           exit iteration_loop
        end if

        if (s(ipt) < t(1)) then
           s(ipt) = t(1)
        end if

        if (s(ipt) > t(nctl+k)) then
           s(ipt) = t(nctl+k)
        end if

        if (norm((s(ipt)-s0(ipt))*deriv,ndim) <= eps1) then
           exit iteration_loop
        end if

        s0(ipt) = s (ipt)
        call eval_curve(s0(ipt),t,k,coef,nctl,ndim,val)
        call eval_curve_deriv(s0(ipt),t,k,coef,nctl,ndim,deriv)
        call eval_curve_deriv2(s0(ipt),t,k,coef,nctl,ndim,deriv2)
        Diff(ipt,:) = val-x0(ipt,:)
        s(ipt) = s0(ipt) - dot_product(deriv,Diff(ipt,:))/(dot_product(deriv2,Diff(ipt,:)) + norm(deriv,ndim)**2)
     end do iteration_loop
  end do
  if (brute_force .eqv. .true.) then
     deallocate(curve_vals)
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
  !     x0      - Real, array size(N,ndim) points we are trying to invert
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     coef    - Real,Array of B-spline coefficients  Size (nctlu,nctlv,ndim)
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
  !     diff    - Real Array size(N,ndim) - Distance between x0 and S(u,v)

  implicit none
  ! Input
  double precision, intent(in)          :: x0(N,ndim)
  integer         , intent(in)          :: ku,kv,nctlu,nctlv,ndim,niter,N
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
  double precision, intent(in)          :: coef(nctlu,nctlv,ndim)
  double precision, intent(in)          :: eps1,eps2
  ! Output
  double precision, intent(inout)       :: u(N),v(N)
  double precision, intent(out)         :: diff(N,ndim)


  ! Working
  double precision                      :: val(ndim),deriv(2,ndim),deriv2(2,2,ndim)
  double precision                      :: val0(ndim)
  logical                               :: brute_force
  integer                               :: i,j,ii,jj,counter,ipt
  double precision                      :: D,D0,u0(N),v0(N),delta(2)
  double precision                      :: A(2,2),ki(2)
  integer                               :: n_sub  ! Huristic Value
  integer                               :: istart,nuu,nvv

  ! Alloctable
  double precision,allocatable          :: surface_vals(:,:,:),uu(:),vv(:)

  ! Functions     
  double precision                      :: norm

  n_sub = 6
  brute_force = .false.
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
     nuu = nctlu-ku+2 + (nctlu-ku+1)*n_sub
     nvv = nctlv-kv+2 + (nctlv-kv+1)*n_sub
     allocate (uu(nuu),vv(nvv))

     counter = 1
     do i=1,nctlu-ku+1
        if (i==1) then
           istart = 0
        else
           istart = 1
        end if
        do j=istart,n_sub+1
           uu(counter) = tu(i+ku-1) + (real(j)/(n_sub+1)*(tu(i+ku)-tu(i+ku-1)))
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
        do j=istart,n_sub+1
           vv(counter) = tv(i+kv-1) + (real(j)/(n_sub+1)*(tv(i+kv)-tv(i+kv-1)))
           counter = counter + 1
        end do
     end do
     allocate(surface_vals(nuu,nvv,ndim))

     do i=1,nuu
        do j=1,nvv
           call eval_surface(uu(i),vv(j),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,&
                surface_vals(i,j,:))
        end do
     end do

     ! Now do comparison
     u0(:) = u(:)   
     v0(:) = v(:)

     do ipt=1,N
        if (u(ipt) < 0 .or. u(ipt) > 1 .or. v(ipt) < 0 .or. v(ipt) > 1) then
           val0 = surface_vals(1,1,:)
           D0 = norm(val0-x0(ipt,:),ndim)
           u0(ipt) = uu(1)
           v0(ipt) = vv(1)
           do i = 1,nuu
              do j = 1,nvv
                 D = norm(surface_vals(i,j,:)-X0(ipt,:),ndim)
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
  do ipt=1,N
     call eval_surface(u0(ipt),v0(ipt),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
     call eval_surface_deriv(u0(ipt),v0(ipt),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv)
     call eval_surface_deriv2(u0(ipt),v0(ipt),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv2)
     Diff(ipt,:) = val-x0(ipt,:)
     u(ipt) = u0(ipt)
     v(ipt) = v0(ipt)
     iteration_loop: do i=1,niter
        ! Check the convergence criteria
        if (norm(Diff(ipt,:),ndim) <= eps1) then
           exit iteration_loop
        end if

        if (norm(dot_product(deriv(1,:),Diff(ipt,:)),ndim)/(norm(deriv(1,:),ndim)*norm(Diff(ipt,:),ndim)) <= eps2 .and. &
             norm(dot_product(deriv(2,:),Diff(ipt,:)),ndim)/(norm(deriv(2,:),ndim)*norm(Diff(ipt,:),ndim)) <= eps2 ) then
           exit iteration_loop
        end if
        u0(ipt) = u(ipt)
        v0(ipt) = v(ipt)
        call eval_surface(u(ipt),v(ipt),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
        call eval_surface_deriv(u(ipt),v(ipt),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv)
        call eval_surface_deriv2(u(ipt),v(ipt),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv2)

        Diff(ipt,:) = val-x0(ipt,:)

        A(1,1) = norm(deriv(1,:),ndim)**2 + dot_product(Diff(ipt,:),deriv2(1,1,:))
        A(1,2) = dot_product(deriv(1,:),deriv(2,:)) + dot_product(Diff(ipt,:),deriv2(1,2,:))
        A(2,1) = A(1,2)
        A(2,2) = norm(deriv(2,:),ndim)**2 + dot_product(Diff(ipt,:),deriv2(2,2,:))

        ki(1) = -dot_product(Diff(ipt,:),deriv(1,:))
        ki(2) = -dot_product(Diff(ipt,:),deriv(2,:))

        call solve_2by2(A,ki,delta)

        u(ipt) = u0(ipt) + delta(1)
        v(ipt) = v0(ipt) + delta(2)

        ! Bounds Checking
        if (u(ipt) < tu(1)) then
           u(ipt) = tu(1)
        end if

        if (u(ipt) > tu(nctlu+ku)) then
           u(ipt) = tu(nctlu+ku)
        end if

        if (v(ipt) < tv(1)) then
           v(ipt) = tv(1)
        end if

        if (v(ipt) > tv(nctlv+kv)) then
           v(ipt) = tv(nctlv+kv)
        end if

        ! No Change convergence Test

        if (norm( (u(ipt)-u0(ipt))*deriv(1,:) + (v(ipt)-v0(ipt))*deriv(2,:),ndim) <= eps1) then
           exit iteration_loop
        end if
     end do iteration_loop
  end do
  if (brute_force .eqv. .true.) then
     deallocate(surface_vals,uu,vv)
  end if
end subroutine point_surface


subroutine point_volume(x0,tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,N,niter,eps1,eps2,u,v,w,Diff)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: pint_surface attempts to solve the point inversion problem for a surface
  !
  !     Description of Arguments
  !     Input
  !     x0      - Real, array size(N,ndim) points we are trying to invert
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     tw      - Real,Knot vector in w. Length nctlw+kw
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     kw      - Integer, order of B-spline in w
  !     coef    - Real,Array of B-spline coefficients  Size (nctlu,nctlv,nctlw,ndim)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     nctlw   - Integer,Number of control points in w
  !     ndim    - Integer, Spatial Dimension
  !     Niter   - Integer, Maximum number of Netwton iterations
  !     eps1    - Real - Eculdian Distance Convergence Measure
  !     eps2    - Real - Cosine Convergence Measure
  !
  !     Ouput 
  !     u       - Real,vector size(N), u parameters where V(u,v,w) is closest to x0
  !     v       - Real,vector size(N), v parameters where V(u,v,w) is closest to x0
  !     w       - Real,vector size(N), w parameters where V(u,v,w) is closest to x0
  !     diff    - Real Array size(N,ndim) - Distance between x0 and V(u,v,w)

  implicit none
  ! Input
  double precision, intent(in)          :: x0(N,ndim)
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,ndim,niter,N
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  double precision, intent(in)          :: coef(nctlu,nctlv,nctlw,ndim)
  double precision, intent(in)          :: eps1,eps2
  ! Output
  double precision, intent(inout)       :: u(N),v(N),w(N)
  double precision, intent(out)         :: diff(N,ndim)


  ! Working
  double precision                      :: val(ndim)
  double precision                      :: deriv(3,ndim)
  double precision                      :: deriv2(3,3,ndim)
  double precision                      :: val0(ndim),uu,vv
  logical                               :: brute_force
  integer                               :: i,j,ii,jj,counter,ipt
  double precision                      :: D,D0,u0(N),v0(N),w0(N),delta(3)
  double precision                      :: A(3,3),ki(3)
  integer                               :: n_sub ! Huristic Value

  ! Alloctable
  double precision,allocatable          :: volume_vals(:,:)

  ! Functions     
  double precision                      :: norm

  n_sub = 3
  brute_force = .false.
  ! First we will evaluate the surface at n points inside each knot span in each direction

  !if we are given a bad guess do the brute force
  point_loop: do ipt=1,N
     if (u(ipt) < 0 .or. u(ipt) > 1 .or. v(ipt) < 0 .or. v(ipt) > 1) then
        brute_force = .True.
        exit point_loop
     end if
  end do point_loop

  if (brute_force .eqv. .True.) then
     print *,'Brute Force Not Implemented Yet'
     stop
  else
     u0(:) = u(:)   
     v0(:) = v(:)
     w0(:) = w(:)
  end if

  ! Now we have u0,v0,w0 so we can do the newton search
  do ipt=1,N
     call eval_volume(u0(ipt),v0(ipt),w0(ipt),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,val)
     call eval_volume_deriv(u0(ipt),v0(ipt),w0(ipt),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,deriv)
     call eval_volume_deriv2(u0(ipt),v0(ipt),w0(ipt),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,deriv2)

     Diff(ipt,:) = val-x0(ipt,:)
     u(ipt) = u0(ipt)
     v(ipt) = v0(ipt)
     w(ipt) = w0(ipt)
     iteration_loop: do i=1,niter
        ! Check the convergence criteria
        if (norm(Diff(ipt,:),ndim) <= eps1) then
           exit iteration_loop
        end if

        if (norm(dot_product(deriv(1,:),Diff(ipt,:)),ndim)/(norm(deriv(1,:),ndim)*norm(Diff(ipt,:),ndim)) <= eps2 .and. &
             norm(dot_product(deriv(2,:),Diff(ipt,:)),ndim)/(norm(deriv(2,:),ndim)*norm(Diff(ipt,:),ndim))<= eps2 .and. & 
             norm(dot_product(deriv(3,:),Diff(ipt,:)),ndim)/(norm(deriv(3,:),ndim)*norm(Diff(ipt,:),ndim))<= eps2) then
           exit iteration_loop
        end if
        u0(ipt) = u(ipt)
        v0(ipt) = v(ipt)
        w0(ipt) = w(ipt)
        call eval_volume(u(ipt),v(ipt),w(ipt),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,val)
        call eval_volume_deriv(u(ipt),v(ipt),w(ipt),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,deriv)
        call eval_volume_deriv2(u(ipt),v(ipt),w(ipt),tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,deriv2)

        Diff(ipt,:) = val-x0(ipt,:)

        A(1,1) = norm(deriv(1,:),ndim)**2 + dot_product(Diff(ipt,:),deriv2(1,1,:))
        A(1,2) = dot_product(deriv(1,:),deriv(2,:)) + dot_product(Diff(ipt,:),deriv2(1,2,:))
        A(1,3) = dot_product(deriv(1,:),deriv(3,:)) + dot_product(Diff(ipt,:),deriv2(1,3,:))
        A(2,2) = norm(deriv(2,:),ndim)**2 + dot_product(Diff(ipt,:),deriv2(2,2,:))
        A(2,3) = dot_product(deriv(2,:),deriv(3,:)) + dot_product(Diff(ipt,:),deriv2(2,3,:))
        A(3,3) = norm(deriv(3,:),ndim)**2 + dot_product(Diff(ipt,:),deriv2(3,3,:))
        A(2,1) = A(1,2)
        A(3,1) = A(1,3)
        A(3,2) = A(2,3)

        ki(1) = -dot_product(Diff(ipt,:),deriv(1,:))
        ki(2) = -dot_product(Diff(ipt,:),deriv(2,:))
        ki(3) = -dot_product(Diff(ipt,:),deriv(3,:))

        call solve_3by3(A,ki,delta)

        u(ipt) = u0(ipt) + delta(1)
        v(ipt) = v0(ipt) + delta(2)
        w(ipt) = w0(ipt) + delta(3)

        ! -- U Bounds Checking --
        if (u(ipt) < tu(1)) then
           u(ipt) = tu(1)
        end if

        if (u(ipt) > tu(nctlu+ku)) then
           u(ipt) = tu(nctlu+ku)
        end if

        ! -- V Bounds Checking --
        if (v(ipt) < tv(1)) then
           v(ipt) = tv(1)
        end if

        if (v(ipt) > tv(nctlv+kv)) then
           v(ipt) = tv(nctlv+kv)
        end if

        ! -- W Bounds Checking --
        if (w(ipt) < tw(1)) then
           w(ipt) = tw(1)
        end if

        if (w(ipt) > tw(nctlw+kw)) then
           w(ipt) = tw(nctlw+kw)
        end if


        ! No Change convergence Test

        if (norm( (u(ipt)-u0(ipt))*deriv(1,:) + &
             (v(ipt)-v0(ipt))*deriv(2,:) + &
             (w(ipt)-w0(ipt))*deriv(3,:),ndim) <= eps1) then
           exit iteration_loop
        end if
     end do iteration_loop
  end do
  if (brute_force .eqv. .true.) then
     deallocate(volume_vals)
  end if
end subroutine point_volume


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
  double precision, intent(in)     :: coef1(n1,ndim),coef2(n2,ndim)
  double precision, intent(inout)  :: s,t 
  integer         , intent(in)     :: Niter
  double precision, intent(in)     :: eps1,eps2

  double precision, intent(out)    :: Diff(ndim)

  ! Working
  double precision                 :: val0(ndim),val1(ndim),val2(ndim)
  double precision                 :: deriv_c1(ndim),deriv2_c1(ndim)
  double precision                 :: deriv_c2(ndim),deriv2_c2(ndim)
  integer                          :: i,j,ii,jj
  double precision                 :: D,D0,s0,t0,delta(2)
  double precision                 :: ki(2),A(2,2)
  integer                          :: n ! Huristic Value

  ! Functions 
  double precision                 :: norm


  ! Is 3 Good Here?? Huristic
  if (k1 == 2 .or. k2 == 2) then
     n = 5
  else
     n = 3
  end if

  if (s < 0 .or. s > 1 .or. t < 0 .or. t > 1) then
     ! Do a brute force approach to get good starting point

     ! Starting point is end of each curve
     call eval_curve(t1(1),t1,k1,coef1,n1,ndim,val1)
     call eval_curve(t2(1),t2,k2,coef2,n2,ndim,val2)
     D0 = norm(val1-val2,ndim)

     do i=1,n1-k1+1
        do ii=1,n
           do j=1,n2-k2+1
              do jj=1,n
                 s = t1(i+k1-1) + (real(ii)/n)*(t1(i+k1)-t1(i+k1-1))
                 t = t2(j+k2-1) + (real(jj)/n)*(t2(j+k2)-t2(j+k2-1))
                 call eval_curve(s,t1,k1,coef1,n1,ndim,val1)
                 call eval_curve(t,t2,k2,coef2,n2,ndim,val2)
                 D = norm(val1-val2,ndim)
                 if (D<D0) then
                    s0 = s
                    t0 = t
                    D0 = D
                 end if
              end do
           end do
        end do
     end do
  else
     s0 = s
     t0 = t
     call eval_curve(s,t1,k1,coef1,n1,ndim,val1)
     call eval_curve(t,t2,k2,coef2,n2,ndim,val2)
     D0 = norm(val1-val2,ndim)
  end if
  ! Now we have s0 and t0 so we can do the newton iteration
  call eval_curve(s,t1,k1,coef1,n1,ndim,val1)
  call eval_curve_deriv(s,t1,k2,coef1,n1,ndim,deriv_c1)
  call eval_curve_deriv2(s,t1,k2,coef1,n1,ndim,deriv2_c1)

  call eval_curve(t,t2,k2,coef2,n2,ndim,val2)
  call eval_curve_deriv(t,t2,k2,coef2,n2,ndim,deriv_c2)
  call eval_curve_deriv2(t,t2,k2,coef2,n2,ndim,deriv2_c2)
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

     call eval_curve(s,t1,k1,coef1,n1,ndim,val1)
     call eval_curve_deriv(s,t1,k2,coef1,n1,ndim,deriv_c1)
     call eval_curve_deriv2(s,t1,k2,coef1,n1,ndim,deriv2_c1)

     call eval_curve(t,t2,k2,coef2,n2,ndim,val2)
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

     s = s0 + delta(1)
     t = t0 + delta(2)

     ! Bounds checking
     if (s < t1(1)) then
        s = t1(1)
     end if

     if (s > t1(n1+k1)) then
        s = t1(n1+k1)
     end if

     if (t < t2(1)) then
        t = t2(1)
     end if

     if (t > t2(n2+k2)) then
        t = t2(n2+k2)
     end if

     ! No change convergence Test
     if (norm( (s-s0)*deriv_c1 + (t-t0)*deriv_c2,ndim) <= eps1) then
        exit iteration_loop
     end if

  end do iteration_loop
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
  !     coefc   - Real, aray size(nctlc,ndim) coefficients for curve
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     coefs   - Real,Array of B-spline coefficients  Size (nctlu,nctlv,ndim)
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
  double precision, intent(in)          :: coefc(nctlc,ndim),coefs(nctlu,nctlv,ndim)
  double precision, intent(in)          :: eps1,eps2
  ! Output
  double precision, intent(inout)       :: u,v,s
  double precision, intent(out)         :: diff(ndim)

  ! Working
  double precision                      :: val_s(ndim),deriv_s(2,ndim),deriv2_s(2,2,ndim)
  double precision                      :: val_c(ndim),deriv_c(ndim),deriv2_C(ndim)
  double precision                      :: val0_s(ndim),val0_c(ndim)
  integer                               :: i,j,ii,jj,l,ll,counter1,counter2
  double precision                      :: D,D0,u0,v0,s0,delta(3)
  double precision                      :: A(3,3),ki(3)
  integer                               :: n,nc ! Huristic Value
  integer                               :: istart,nss,nuu,nvv
  logical                               :: brute_force
  ! Allocatable 
  double precision,allocatable          :: curve_vals(:,:),uu(:),vv(:),ss(:)
  double precision,allocatable          :: surface_vals(:,:,:)

  ! Functions     
  double precision                      :: norm

  n = 6
  nc = 12
  ! First we will evaluate the surface at n points inside each knot span in each direction
  !if we are given a bad guess do the brute force

  if (u < 0 .or. u > 1 .or. v < 0 .or. v > 1 .or. s < 0 .or. s > 1) then
     brute_force = .True.

     ! Generate the ss,uu and vv values
     nss = nctlc-kc+2 + (nctlc-kc+1)*nc
     nuu = nctlu-ku+2 + (nctlu-ku+1)*n
     nvv = nctlv-kv+2 + (nctlv-kv+1)*n
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
        do j=istart,n+1
           uu(counter1) = tu(i+ku-1) + (real(j)/(n+1)*(tu(i+ku)-tu(i+ku-1)))
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
        do j=istart,n+1
           vv(counter1) = tv(i+kv-1) + (real(j)/(n+1)*(tv(i+kv)-tv(i+kv-1)))
           counter1 = counter1 + 1
        end do
     end do

     allocate(curve_vals(nss,ndim),surface_vals(nuu,nvv,ndim))

     do i=1,nss
        call eval_curve(ss(i),tc,kc,coefc,nctlc,ndim,curve_vals(i,:))
     end  do

     do i=1,nuu
        do j=1,nvv
           call eval_surface(uu(i),vv(j),tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,&
                surface_vals(i,j,:))
        end do
     end do

     ! Now do comparison
     val0_s = surface_vals(1,1,:)
     val0_c = curve_vals(1,:)
     D0 = norm(val0_s-val0_c,ndim)
     u0 = uu(1)
     v0 = vv(1)
     s0 = ss(1)
     do l = 1,nss
        do i = 1,nuu
           do j = 1,nvv
              D = norm(surface_vals(i,j,:)-curve_vals(l,:),ndim)
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

  call eval_surface(u0,v0,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,val_s)
  call eval_surface_deriv(u0,v0,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,deriv_s)
  call eval_surface_deriv2(u0,v0,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,deriv2_s)

  call eval_curve(s0,tc,kc,coefc,nctlc,ndim,val_c)
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

     if ( norm(dot_product(deriv_s(1,:),Diff),ndim)/(norm(deriv_s(1,:),ndim)*norm(Diff,ndim)) <= eps2 .and. &
          norm(dot_product(deriv_s(2,:),Diff),ndim)/(norm(deriv_s(2,:),ndim)*norm(Diff,ndim)) <= eps2 .and. &
          norm(dot_product(deriv_c     ,Diff),ndim)/(norm(deriv_c     ,ndim)*norm(Diff,ndim)) <= eps2) then
        exit iteration_loop
     end if
     u0 = u
     v0 = v
     s0 = s

     call eval_surface(u0,v0,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,val_s)
     call eval_surface_deriv(u0,v0,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,deriv_s)
     call eval_surface_deriv2(u0,v0,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,deriv2_s)

     call eval_curve(s0,tc,kc,coefc,nctlc,ndim,val_c)
     call eval_curve_deriv(s0,tc,kc,coefc,nctlc,ndim,deriv_c)
     call eval_curve_deriv2(s0,tc,kc,coefc,nctlc,ndim,deriv2_c)

     Diff = val_s-val_c

     A(1,1) = norm(deriv_s(1,:),ndim)**2 + dot_product(Diff,deriv2_s(1,1,:))
     A(1,2) = dot_product(deriv_s(1,:),deriv_s(2,:)) + dot_product(Diff,deriv2_s(1,2,:))
     A(1,3) = -dot_product(deriv_c,deriv_s(1,:))

     A(2,1) = A(1,2)
     A(2,2) = norm(deriv_s(2,:),ndim)**2 + dot_product(Diff,deriv2_s(2,2,:))
     A(2,3) = -dot_product(deriv_c,deriv_s(2,:))

     A(3,1) = dot_product(deriv_s(1,:),deriv_c) + dot_product(Diff,deriv2_c)
     A(3,2) = dot_product(deriv_s(2,:),deriv_c) + dot_product(Diff,deriv2_c)
     A(3,3) = -norm(deriv_c,ndim)**2 + dot_product(Diff,deriv2_c)

     ki(1) = -dot_product(Diff,deriv_s(1,:))
     ki(2) = -dot_product(Diff,deriv_s(2,:))
     ki(3) = -dot_product(Diff,deriv_c)
     call solve_3by3(A,ki,delta)

     u = u0 + delta(1)
     v = v0 + delta(2)
     s = s0 + delta(3)
     ! Bounds Checking
     if (u < tu(1)) then
        u = tu(1)
     end if

     if (u > tu(nctlu+ku)) then
        u = tu(nctlu+ku)
     end if

     if (v < tv(1)) then
        v = tv(1)
     end if

     if (v > tv(nctlv+kv)) then
        v = tv(nctlv+kv)
     end if

     if (s < tc(1)) then
        s = tc(1)
     end if

     if (s > tc(nctlc+kc)) then
        s = tc(nctlc+kc)
     end if

     ! No Change convergence Test

     if (norm((u-u0)*deriv_s(1,:)+(v-v0)*deriv_s(2,:)+(s-s0)*deriv_c,ndim) <= eps1) then
        exit iteration_loop
     end if

  end do iteration_loop

  if (brute_force .eqv. .True.) then
     deallocate(curve_vals,surface_vals,ss,uu,vv)
  end if

end subroutine curve_surface

subroutine solve_2by2(A,b,x)
  ! Solve a 2 x 2 system  -- With NO checking
  double precision         :: A(2,2),b(2),x(2)
  double precision         :: idet

  idet = 1/(A(1,1)*A(2,2)-A(1,2)*A(2,1))
  x(1) = idet*(A(2,2)*b(1) - A(1,2)*b(2))
  x(2) = idet*(-A(2,1)*b(1) + A(1,1)*b(2))

end subroutine solve_2by2

subroutine solve_3by3(A,b,x)
  ! Solve a 3 x 3 system  -- With NO checking
  double precision         :: A(3,3),b(3),x(3)
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
