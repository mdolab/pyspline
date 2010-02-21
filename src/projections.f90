subroutine point_curve(x0,t,k,coef,nctl,ndim,Niter,eps1,eps2,s,Diff)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: point_curve attempts to solve the point inversion problem
  !
  !     Description of Arguments
  !     Input
  !     x0      - Real, array size(ndim) point we are trying to invert
  !     s0      - Real, initial guess for the solution
  !     t       - Real, Knot vector. Length nctl+k
  !     k       - Integer,order of B-spline
  !     coef    - Real,Array of B-spline coefficients and weights. Size (nctl,ndim)
  !     nctl    - Integer,Number of control points
  !     ndim    - Integer, spatial dimension of curve
  !     Niter   - Integer, Maximum number of Netwton iterations
  !     eps1    - Real - Eculdian Distance Convergence Measure
  !     eps2    - Real - Cosine Convergence Measure
  !     s       - Real, guess parameter where C(s) is closest to x0 
  !
  !     Ouput 
  !     s       - Real, parameter where C(s) is closest to x0
  !     diff    - Real Array size(ndim) - Distance between x0 and curve(s)

  implicit none
  ! Input
  double precision, intent(in)          :: x0(ndim)
  integer         , intent(in)          :: k,nctl,ndim
  double precision, intent(in)          :: t(nctl+k)
  double precision, intent(in)          :: coef(nctl,ndim)
  double precision, intent(in)          :: Niter
  double precision, intent(in)          :: eps1,eps2

  ! Output
  double precision, intent(out)         :: s,Diff(ndim)

  ! Working
  double precision                      :: val(ndim),deriv(ndim),deriv2(ndim)
  double precision                      :: val0(ndim),s0
  integer                               :: i,j
  double precision                      :: D,D0
  integer                               :: n ! Huristic Value

  ! Functions
  double precision                      :: norm

  ! Is 3  Good Here?? Huristic!
  n=10

  !if we are given a bad guess do the brute force
  if (s < 0 .or. s > 1) then
     ! First we evaluate the curve at n points inside each knot span
     call eval_curve(t(1),t,k,coef,nctl,ndim,val0)
     D0 = norm(val0-x0,ndim)
     do i=1,nctl-k+1 ! Number of knot spans for non-a-knot end conditions
        do j = 1,n
           s = t(i+k-1) + (real(j)/n)*(t(i+k)-t(i+k-1))
           call eval_curve(s,t,k,coef,nctl,ndim,val)
           D  = norm(x0-val,ndim)
           if (D<D0) then
              s0 = s
              D0 = D
           end if
        end do
     end do
  else
     s0 = s
     call eval_curve(s0,t,k,coef,nctl,ndim,val)
     D0 = norm(x0-val,ndim)
  end if

  ! Now we have s0 so we should be able to do the netwon seach
  call eval_curve(s0,t,k,coef,nctl,ndim,val)
  call eval_curve_deriv(s0,t,k,coef,nctl,ndim,deriv)
  call eval_curve_deriv2(s0,t,k,coef,nctl,ndim,deriv2)
  Diff = val-x0
  s = s0 - dot_product(deriv,Diff)/(dot_product(deriv2,Diff) + norm(deriv,ndim)**2)
  iteration_loop: do i =1,Niter
     !print *,'iter',i,s0,s
     ! Check the Convergence Criteria
     if (norm(Diff,ndim) <= eps1) then
        exit iteration_loop
     end if

     if (  norm(dot_product(deriv,Diff),ndim)/ (norm(deriv,ndim)*norm(Diff,ndim)) <= eps2) then
        exit iteration_loop
     end if

     if (s < t(1)) then
        s = t(1)
     end if

     if (s > t(nctl+k)) then
        s = t(nctl+k)
     end if

     if (norm((s-s0)*deriv,ndim) <= eps1) then
        exit iteration_loop
     end if

     s0 = s 
     call eval_curve(s0,t,k,coef,nctl,ndim,val)
     call eval_curve_deriv(s0,t,k,coef,nctl,ndim,deriv)
     call eval_curve_deriv2(s0,t,k,coef,nctl,ndim,deriv2)
     Diff = val-x0
     s = s0 - dot_product(deriv,Diff)/(dot_product(deriv2,Diff) + norm(deriv,ndim)**2)

  end do iteration_loop

end subroutine point_curve

subroutine point_surface(x0,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,niter,eps1,eps2,u,v,Diff)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: pint_surface attempts to solve the point inversion problem for a surface
  !
  !     Description of Arguments
  !     Input
  !     x0      - Real, array size(ndim) point we are trying to invert
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
  !     u       - Real, u parameter where S(u,v) is closest to x0
  !     v       - Real, v parameter where S(u,v) is closest to x0
  !     diff    - Real Array size(ndim) - Distance between x0 and S(u,v)

  implicit none
  ! Input
  double precision, intent(in)          :: x0(ndim)
  integer         , intent(in)          :: ku,kv,nctlu,nctlv,ndim,niter
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
  double precision, intent(in)          :: coef(nctlu,nctlv,ndim)
  double precision, intent(in)          :: eps1,eps2
  ! Output
  double precision, intent(inout)       :: u,v
  double precision, intent(out)         :: diff(ndim)


  ! Working
  double precision                      :: val(ndim),deriv(2,ndim),deriv2(2,2,ndim)
  double precision                      :: val0(ndim)
  integer                               :: i,j,ii,jj
  double precision                      :: D,D0,u0,v0,delta(2)
  double precision                      :: A(2,2),ki(2)
  integer                               :: n ! Huristic Value

  ! Functions     
  double precision                      :: norm

  n = 3

  ! First we will evaluate the surface at n points inside each knot span in each direction

  !if we are given a bad guess do the brute force
  if (u < 0 .or. u > 1 .or. v < 0 .or. v > 1) then


     call eval_surface(tu(1),tv(1),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val0)
     D0 = norm(val0-x0,ndim)

     do i = 1,nctlu-ku+1
        do ii = 1,n
           do j = 1,nctlv-kv+1
              do jj =1,n
                 u = tu(i+ku-1) + (real(ii)/n)*(tu(i+ku)-tu(i+ku-1))
                 v = tv(j+kv-1) + (real(jj)/n)*(tv(i+kv)-tv(i+kv-1))
                 call eval_surface(u,v,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
                 D = norm(x0-val,ndim)
                 if (D<D0) then
                    u0 = u
                    v0 = v
                    D0 = D
                 end if
              end do
           end do
        end do
     end do
  else
     u0 = u
     v0 = v
     call eval_surface(u,v,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
     D0 = norm(x0-val,ndim)
  end if
  ! Now we have u0 and v0 so we can do the newton 

  call eval_surface(u0,v0,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
  call eval_surface_deriv(u0,v0,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv)
  call eval_surface_deriv2(u0,v0,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv2)
  Diff = val-x0

  u = u0
  v = v0
  iteration_loop: do i=1,niter
     !print *,i,u,v
     ! Check the convergence criteria
     if (norm(Diff,ndim) <= eps1) then
        exit iteration_loop
     end if

     if (norm(dot_product(deriv(1,:),Diff),ndim)/(norm(deriv(1,:),ndim)*norm(Diff,ndim)) <= eps2 .and. &
          norm(dot_product(deriv(2,:),Diff),ndim)/(norm(deriv(2,:),ndim)*norm(Diff,ndim)) <= eps2 ) then
        exit iteration_loop
     end if
     u0 = u
     v0 = v
     call eval_surface(u,v,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
     call eval_surface_deriv(u,v,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv)
     call eval_surface_deriv2(u,v,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv2)

     Diff = val-x0

     A(1,1) = norm(deriv(1,:),ndim)**2 + dot_product(Diff,deriv2(1,1,:))
     A(1,2) = dot_product(deriv(1,:),deriv(2,:)) + dot_product(Diff,deriv2(1,2,:))
     A(2,1) = A(1,2)
     A(2,2) = norm(deriv(2,:),ndim)**2 + dot_product(Diff,deriv2(2,2,:))

     ki(1) = -dot_product(Diff,deriv(1,:))
     ki(2) = -dot_product(Diff,deriv(2,:))

     call solve_2by2(A,ki,delta)

     u = u0 + delta(1)
     v = v0 + delta(2)

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

     if (v > tu(nctlv+kv)) then
        v = tv(nctlv+kv)
     end if

     ! No Change convergence Test

     if (norm( (u-u0)*deriv(1,:) + (v-v0)*deriv(2,:),ndim) <= eps1) then
        exit iteration_loop
     end if

  end do iteration_loop
end subroutine point_surface


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
  integer                               :: i,j,ii,jj,l,ll
  double precision                      :: D,D0,u0,v0,s0,delta(3)
  double precision                      :: A(3,3),ki(3)
  integer                               :: n ! Huristic Value

  ! Functions     
  double precision                      :: norm

  n = 3

  ! First we will evaluate the surface at n points inside each knot span in each direction
  !if we are given a bad guess do the brute force
  if (u < 0 .or. u > 1 .or. v < 0 .or. v > 1 .or. s < 0 .or. s > 1) then
     call eval_surface(tu(1),tv(1),tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,val0_s)
     call eval_curve(tc(1),tc,kc,coefc,nctlc,ndim,val0_c)
     D0 = norm(val0_s-val0_c,ndim)

     do l = 1,nctlc-kc+1
        do ll = 1,n
           do i = 1,nctlu-ku+1
              do ii = 1,n
                 do j = 1,nctlv-kv+1
                    do jj =1,n
                       u = tu(i+ku-1) + (real(ii)/n)*(tu(i+ku)-tu(i+ku-1))
                       v = tv(j+kv-1) + (real(jj)/n)*(tv(i+kv)-tv(i+kv-1))
                       s = tc(l+kc-1) + (real(ll)/n)*(tc(l+kc)-tc(i+kc-1))
                       call eval_surface(u,v,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,val_s)
                       call eval_curve(s,tc,kc,coefc,nctlc,ndim,val_c)
                       D = norm(val_s-val_c,ndim)
                       if (D<D0) then
                          u0 = u
                          v0 = v
                          s0 = s
                          D0 = D
                       end if
                    end do
                 end do
              end do
           end do
        end do
     end do
  else
     u0 = u
     v0 = v
     s0 = s
     call eval_surface(u,v,tu,tv,ku,kv,coefs,nctlu,nctlv,ndim,val_s)
     call eval_curve(s,tc,kc,coefc,nctlc,ndim,val_c)
     D0 = norm(val_s-val_c,ndim)
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
     !print *,'u,v,s:',u,v,s,Diff
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

     if (v > tu(nctlv+kv)) then
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
