! This file contains all the projection functionality used in
! pySpline. There are 5 combinations that can yield single solutions. They are:
! 1. point-curve   (1 dof)
! 2. point-surface (2 dof)
! 3. point-volume  (3 dof)
! 4. curve-cruve   (2 dof)
! 5. curve-surface (3 dof)

! Additionally, each combination requires a globlalization function to
! obtain a good startin point for the Newton search if one is not
! already provided. The five combinations of these brute-force
! starting points are also included. 

! Define some parameters used for all projection functions 
#define LSFailMax 2
#define wolfe .001
#define nLine 20

subroutine point_curve(x0, t, k, coef, nctl, ndim, Niter, eps, s, Diff)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: point_curve attempts to solve the point inversion problem
  !
  !     Description of Arguments
  !     Input
  !     x0      - Real, array size(ndim) point we are trying to invert
  !     t       - Real, Knot vector. Length nctl+k
  !     k       - Integer, order of B-spline
  !     coef    - Real, Array of B-spline coefficients and weights. Size (ndim, nctl)
  !     nctl    - Integer, Number of control points
  !     ndim    - Integer, spatial dimension of curve
  !     Niter   - Integer, Maximum number of Netwton iterations
  !     eps     - Real - Eculdian Distance Convergence Measure
  !     s       - Real, vector, length(N), guess parameters where C(s)
  !               is closest to x0 
  !
  !     Ouput 
  !     s       - Real, parameter where C(s) is closest to x0
  !     diff    - Real, array, size(ndim)- Distance between x0 and curve(s)

  use precision
  implicit none

  ! Input
  integer            , intent(in)          :: k, nctl, ndim, niter
  real(kind=realType), intent(in)          :: x0(ndim)
  real(kind=realType), intent(in)          :: t(nctl+k)
  real(kind=realType), intent(in)          :: coef(ndim, nctl)
  real(kind=realType), intent(in)          :: eps

  ! Output
  real(kind=realType), intent(out)         :: s, Diff(ndim)

  ! Working
  real(kind=realType)   :: val(ndim), deriv(ndim), deriv2(ndim), step, c, dist, p_diff
  real(kind=realType)   :: grad, hessian,  update, R(ndim), nDist, fval, nfval, pgrad, newPt
  integer               :: m, ii, NLSFail
  logical               :: flag, cflag

  NLSFail = 0
  iteration_loop: do ii=1, Niter

     ! Evaluate point and its derivatives
     call eval_curve(s, t, k, coef, nctl, ndim, 1, val)
     call eval_curve_deriv(s, t, k, coef, nctl, ndim, deriv)
     call eval_curve_deriv2(s, t, k, coef, nctl, ndim, deriv2)
     
     ! Distance is R, "function value" fval is what we minimize
     R = val - x0
     nDist = NORM2(R)
     fval = 0.5*nDist**2

     ! Calculate the Gradient
     grad = dot_product(R, deriv)

     ! Calculate the Hessian
     hessian = dot_product(deriv,deriv) + dot_product(R, deriv2)

     ! Bounds checking
     flag = .False.
     if (s < t(1) + eps .and. grad >= 0.0) then
        flag = .True. 
        s = t(1)
     end if

     if (s > t(nctl + k)-eps .and. grad <= 0.0) then
        flag = .True. 
        s = t(nctl + k)
     end if

     if (flag) then
        grad = 0.0
        hessian = 1.0
     end if

     ! Check the norm of the gradient
     if (abs(grad) <  eps) then
        exit iteration_loop
     end if

     ! "Invert" the hessian 
     update = grad/hessian
     update = -update
     pgrad = update*grad !dot_product(update, grad)

     ! Check that this is the descent direction
     if (pgrad >= 0.0) then
        update = -grad/ABS(grad)
        pgrad = update*grad !dot_product(update, grad)
     end if

     step = 1.0
     nDist = 0.0
     lineLoop: do m=1, nLine
        newPt = s + step*update
        cflag = .False. ! Check if the constraint is applied

        if (newpt > t(nctl+k)) then
           cflag = .True. 
           newPt = t(nctl+k)
        end if
        if (newPt < t(1)) then
           cflag = .True. 
           newpt = t(1)
        end if
        
        ! Evaluate the new point
        call eval_curve(newPt, t, k, coef, nctl, ndim, 1, val)

        ! Distance is R, "function value" fval is what we minimize
        R = val - x0
        nDist = NORM2(R)
        nfVal = 0.5*nDist**2
        
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

     if ( m == nLine+1 ) then
        dist = ndist
        nLSFail = nLSFail + 1

        if (NLSFail > LSFailMax) then ! There is nothing more we can do...
           exit Iteration_loop
        end if
     else
        NLSFail = 0
        ! Check if there has been no change in the coordinates
        p_diff = ABS(s-newpt)
        if (p_diff < eps) then
           exit Iteration_loop
        end if
     end if

     s = newpt
  end do iteration_loop

  diff = R

end subroutine point_curve

subroutine point_surface(x0, tu, tv, ku, kv, coef, nctlu, nctlv, ndim, niter, eps, u, v, Diff)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: point_surface attempts to solve the point inversion problem for a surface
  !
  !     Description of Arguments
  !     Input
  !     x0      - Real, array size(ndim) point we are trying to invert
  !     tu      - Real, Knot vector in u. Length nctlu+ku
  !     tv      - Real, Knot vector in v. Length nctlv+kv
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     coef    - Real, Array of B-spline coefficients  Size (ndim, nctlv, nctlu)
  !     nctlu   - Integer, Number of control points in u
  !     nctlv   - Integer, Number of control points in v
  !     ndim    - Integer, Spatial Dimension
  !     Niter   - Integer, Maximum number of Netwton iterations
  !     eps     - Real - Eculdian Distance Convergence Measure
  !
  !     Ouput 
  !     u       - Real, u parameter where S(u, v) is closest to x0
  !     v       - Real, v parameter where S(u, v) is closest to x0
  !     diff    - Real Array size(ndim) - Distance between x0 and S(u, v)

  use precision
  implicit none

  ! Input
  integer            , intent(in)     :: ku, kv, nctlu, nctlv, ndim, niter
  real(kind=realType), intent(in)     :: x0(ndim)
  real(kind=realType), intent(in)     :: tu(nctlu+ku), tv(nctlv+kv)
  real(kind=realType), intent(in)     :: coef(ndim, nctlv, nctlu)
  real(kind=realType), intent(in)     :: eps

  ! Output
  real(kind=realType), intent(inout)  :: u, v
  real(kind=realType), intent(out)    :: diff(ndim)

  ! Working
  real(kind=realType)   :: val(ndim), deriv(ndim, 2), deriv2(ndim, 2, 2)
  real(kind=realType)   :: R(nDim), low(2), high(2), pt(2), newpt(2), update(2)
  real(kind=realType)   :: step, dist, nDist, pgrad, fval, nfval, c, p_diff
  real(kind=realType)   :: grad(2), hessian(2, 2)
  integer               :: i, j, m, NLSFail, ii
  logical               :: flag, cflag

  NLSFail = 0

  ! Set lower and upper bounds for u, v based on knot vector
  low(1) = tu(1)
  low(2) = tv(1)
  high(1) = tu(Nctlu+ku)
  high(2) = tv(Nctlv+kv)

  pt(1) = u
  pt(2) = v

  iteration_loop: do ii=1, niter
     call eval_surface       (pt(1), pt(2), tu, tv, ku, kv, coef, nctlu, nctlv, ndim, 1, 1, val)
     call eval_surface_deriv (pt(1), pt(2), tu, tv, ku, kv, coef, nctlu, nctlv, ndim, deriv)
     call eval_surface_deriv2(pt(1), pt(2), tu, tv, ku, kv, coef, nctlu, nctlv, ndim, deriv2)
        
     ! Distance is R, "function value" fval is what we minimize
     R = val - X0
     nDist = NORM2(R)
     fval = 0.5*nDist**2
     
     ! Calculate the Gradient
     do i=1,2
        grad(i) = dot_product(R, deriv(:, i))
     end do
     
     ! Calculate the Hessian
     do j=1, 2
        do i=1, 2
           hessian(i, j) = dot_product(deriv(:, i), deriv(:, j)) + &
                dot_product(R, deriv2(:, i, j))
        end do
     end do

     ! Bounds Checking
     do i=1, 2
        flag = .False.
        if (pt(i) < low(i)+eps .and. grad(i) >= 0.0) then
           flag = .True.
           pt(i) = low(i)
        end if

        if (pt(i) > high(i)-eps .and. grad(i) <= 0.0) then
           flag = .True.
           pt(i) = high(i)
        end if

        if ( flag ) then
           grad(i) = 0.0
           hessian(:, i) = 0.0
           hessian(i, :) = 0.0
           hessian(i, i) = 1.0
        end if
     end do

     if ( NORM2(grad) < eps) then
        exit iteration_loop
     end if
        
     ! Invert the hessian, compute the update and the projected gradient
     call solve_2by2(hessian, grad, update)
     update = -update
     pgrad = dot_product(update, grad)

     !Check that this is a descent direction - 
     !otherwise use the negative gradient    
     if ( pgrad >= 0.0 ) then
        update = -grad/NORM2(grad)
        pgrad = dot_product(update, grad)
     end if
        
     step = 1.0
     nDist = 0.0
     lineloop: do m=1, nLine
          
        newpt = pt + step * update

        cflag = .False. ! Check if the constraints are applied
        ! Check if the new point exceeds the bounds
        do i=1, 2
           if (newpt(i) > high(i)) then
              cflag = .True.
              newpt(i) = high(i)
           else if (newpt(i) < low(i)) then
              cflag = .True.
              newpt(i) = low(i)
           end if
        end do
           
        ! Re-evaluate new position
        call eval_surface(newpt(1), newpt(2), tu, tv, ku, kv, coef, nctlu, nctlv, ndim, 1, 1, val)

        ! Check if the new function value is lower, 
        ! otherwise adjust the step size
        R = val - X0
        ndist = NORM2(R)
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

     if ( m == nLine + 1 ) then
        dist = ndist
        nLSFail = nLSFail + 1
        if (NLSFail > LSFailMax) then 
           exit Iteration_loop
        end if
     else
        NLSFail = 0
        ! Check if there has been no change in the coordinates
        p_diff = NORM2(pt-newpt)
        if (p_diff < eps) then
           exit Iteration_loop
        end if
     end if

     pt = newpt
  end do iteration_loop
     
  ! Set the final values of the parameters and our distance function
  u = pt(1)  
  v = pt(2)

  ! Get the ACTUAL difference
  call eval_surface(pt(1), pt(2), tu, tv, ku, kv, coef, nctlu, nctlv, &
       ndim, 1, 1, val)
  diff = val - X0
    
end subroutine point_surface

subroutine point_volume(x0, tu, tv, tw, ku, kv, kw, coef, nctlu, nctlv, nctlw, ndim, &
     niter, eps, u, v, w, Diff)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: point_volume attempts to solve the point inversion problem for a volume 
  !
  !     Description of Arguments
  !     Input
  !     x0      - Real, array size(ndim) point we are trying to invert
  !     tu      - Real, Knot vector in u. Length nctlu+ku
  !     tv      - Real, Knot vector in v. Length nctlv+kv
  !     tw      - Real, Knot vector in w. Length nctlw+kw
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     kw      - Integer, order of B-spline in w
  !     coef    - Real, Array of B-spline coefficients  Size (ndim, nctlw, nctlv, nctlu)
  !     nctlu   - Integer, Number of control points in u
  !     nctlv   - Integer, Number of control points in v
  !     nctlw   - Integer, Number of control points in w
  !     ndim    - Integer, Spatial Dimension
  !     Niter   - Integer, Maximum number of Netwton iterations
  !     eps     - Real - Eculdian Distance Convergence Measure
  !
  !     Ouput 
  !     u       - Real, u parameters where V(u, v, w) is closest to x0
  !     v       - Real, v parameters where V(u, v, w) is closest to x0
  !     w       - Real, w parameters where V(u, v, w) is closest to x0
  !     diff    - Real Array size(N, ndim) - Distance between x0 and V(u, v, w)

  use precision
  implicit none
  ! Input
  integer            , intent(in)          :: ku, kv, kw, nctlu, nctlv, nctlw, ndim, niter
  real(kind=realType), intent(in)          :: x0(ndim)
  real(kind=realType), intent(in)          :: tu(nctlu+ku), tv(nctlv+kv), tw(nctlw+kw)
  real(kind=realType), intent(in)          :: coef(ndim, nctlw, nctlv, nctlu)
  real(kind=realType), intent(in)          :: eps

  ! Output
  real(kind=realType), intent(inout)       :: u, v, w
  real(kind=realType), intent(out)         :: diff(ndim)

  ! Working
  real(kind=realType)   :: val(ndim), deriv(ndim, 3), deriv2(ndim, 3, 3)
  real(kind=realType)   :: R(nDim), low(3), high(3), pt(3), newpt(3), update(3)
  real(kind=realType)   :: step, dist, nDist, pgrad, fval, nfval, c, p_diff
  real(kind=realType)   :: grad(3), hessian(3, 3)
  integer               :: i, j, m, NLSFail, ii
  logical               :: flag, cflag

  NLSFail = 0
  ! Set lower and upper bounds for u, v, w based on knot vector
  low(1) = tu(1)
  low(2) = tv(1)
  low(3) = tw(1)
  high(1) = tu(Nctlu+ku)
  high(2) = tv(Nctlv+kv)
  high(3) = tw(Nctlw+kw)

  pt(1) = u
  pt(2) = v
  pt(3) = w
  
  iteration_loop: do ii=1, niter
     call eval_volume(pt(1), pt(2), pt(3), tu, tv, tw, ku, kv, kw, coef, &
          nctlu, nctlv, nctlw, ndim, 1, 1, 1, val)
     call eval_volume_deriv(pt(1), pt(2), pt(3), tu, tv, tw, ku, kv, kw, coef, &
          nctlu, nctlv, nctlw, ndim, deriv)
     call eval_volume_deriv2(pt(1), pt(2), pt(3), tu, tv, tw, ku, kv, kw, coef, &
          nctlu, nctlv, nctlw, ndim, deriv2)
        
     ! Distance is R, "function value" fval is what we minimize
     R = val - X0
     nDist = NORM2(R)
     fval = 0.5*nDist**2
     ! Calculate the Gradient
     do i=1,3
        grad(i) = dot_product(R, deriv(:, i))
     end do

     ! Calculate the Hessian
     do j=1, 3
        do i=1, 3
           hessian(i, j) = dot_product(deriv(:, i), deriv(:, j)) + &
                dot_product(R, deriv2(:, i, j))
        end do
     end do
     
     ! Bounds Checking
     do i=1, 3
        flag = .False.
        if (pt(i) < low(i)+eps .and. grad(i) >= 0.0) then
           flag = .True.
           pt(i) = low(i)
        end if
        
        if (pt(i) > high(i)-eps .and. grad(i) <= 0.0) then
           flag = .True.
           pt(i) = high(i)
        end if

        if ( flag ) then
           grad(i) = 0.0
           hessian(:, i) = 0.0
           hessian(i, :) = 0.0
           hessian(i, i) = 1.0
        end if
     end do

     if (NORM2(grad) < eps) then
        exit iteration_loop
     end if
        
     ! Invert the hessian, compute the update and the projected gradient
     call solve_3by3(hessian, grad, update)

     update = -update
     pgrad = dot_product(update, grad)
     
     !Check that this is a descent direction - 
     !otherwise use the negative gradient    
     if ( pgrad >= 0.0 ) then
        update = -grad/NORM2(grad)
        pgrad = dot_product(update, grad)
     end if
        
     step = 1.0
     nDist = 0.0
     lineloop: do m=1, nLine
          
        newpt = pt + step * update
           
        cflag = .False. ! Check if the constraints are applied
        ! Check if the new point exceeds the bounds
        do i=1, 3
           if (newpt(i) > high(i)) then
              cflag = .True.
              newpt(i) = high(i)
           else if (newpt(i) < low(i)) then
              cflag = .True.
              newpt(i) = low(i)
           end if
        end do
           
        ! Re-evaluate new position
        call eval_volume(newpt(1), newpt(2), newpt(3), &
             tu, tv, tw, ku, kv, kw, coef, nctlu, nctlv, nctlw, ndim, 1, 1, 1, val)
        
        ! Check if the new function value is lower, 
        ! otherwise adjust the step size
        R = val - X0
        ndist = NORM2(R)
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

     if ( m == nLine + 1 ) then
        dist = ndist
        nLSFail = nLSFail + 1
        if (NLSFail > LSFailMax) then ! There is nothing more we can do...
           exit Iteration_loop
        end if
     else
        NLSFail = 0
        ! Check if there has been no change in the coordinates
        p_diff = NORM2(pt-newpt)
        if (p_diff < eps) then
           exit Iteration_loop
        end if
     end if
     pt = newpt

  end do iteration_loop
     
  ! Set the final values of the parameters and our distance function
  u = pt(1)  
  v = pt(2)
  w = pt(3)

  ! Get the ACTUAL difference
  call eval_volume(pt(1), pt(2), pt(3), &
       tu, tv, tw, ku, kv, kw, coef, nctlu, nctlv, nctlw, ndim, 1, 1, 1, val)
  diff = val - X0
    
end subroutine point_volume

subroutine curve_curve(t1, k1, coef1, t2, k2, coef2, n1, n2, ndim, Niter,&
     eps, s, t, Diff)

  !*** DESCRIPTION
  !
  !     Written by Gaetan Kenway
  ! 
  !     The function takes the spline defination of a bivariate spline
  !     surface as defined by coef, kx, ky, nx, ny, tx, ty and a point x0 and
  !     determines the u, v positions that minimizes the distance between
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
  !     Diff    - Real: Distance between curves
  !

  use precision
  implicit none

  ! Input/Output
  integer            , intent(in)     :: n1, n2, k1, k2, ndim
  real(kind=realType), intent(in)     :: t1(n1+k1), t2(n2+k2)
  real(kind=realType), intent(in)     :: coef1(ndim, n1), coef2(ndim, n2)
  real(kind=realType), intent(inout)  :: s, t 
  integer            , intent(in)     :: Niter
  real(kind=realType), intent(in)     :: eps
  real(kind=realType), intent(out)    :: Diff(ndim)

  ! Working
  real(kind=realType) :: val1(ndim), val2(ndim)
  real(kind=realType) :: deriv_c1(ndim), deriv2_c1(ndim)
  real(kind=realType) :: deriv_c2(ndim), deriv2_c2(ndim)
  integer             :: i, m, ii, NLSFail
  real(kind=realType) :: low(2), high(2), pt(2), R(ndim), ndist
  real(kind=realType) :: fval, nfval, dist
  real(kind=realType) :: hessian(2, 2), grad(2), newpt(2), c
  real(kind=realType) :: pgrad, update(2), step, p_diff
  logical :: flag, cflag
  
  NLSFail = 0

  low(1) = t1(1)
  low(2) = t2(1)
  high(1) = t1(n1+k1)
  high(2) = t2(n2+k2)

  pt(1) = s
  pt(2) = t

  iteration_loop: do ii=1, niter

     call eval_curve       (pt(1), t1, k1, coef1, n1, ndim, 1, val1)
     call eval_curve_deriv (pt(1), t1, k2, coef1, n1, ndim, deriv_c1)
     call eval_curve_deriv2(pt(1), t1, k2, coef1, n1, ndim, deriv2_c1)
     
     call eval_curve       (pt(2), t2, k2, coef2, n2, ndim, 1, val2)
     call eval_curve_deriv (pt(2), t2, k2, coef2, n2, ndim, deriv_c2)
     call eval_curve_deriv2(pt(2), t2, k2, coef2, n2, ndim, deriv2_c2)

     ! Distance is R, "function value" fval is what we minimize
     R = val1-val2
     nDist = NORM2(R)
     fval = 0.5*nDist**2

     ! Calculate the Gradient
     grad(1) = dot_product(R, deriv_c1)
     grad(2) = dot_product(R, -deriv_c2)

     ! Calculate the Hessian
     hessian(1, 1) = dot_product(deriv_c1, deriv_c1) + dot_product(R, deriv2_c1)
     hessian(1, 2) = dot_product(deriv_c1, -deriv_c2)
     hessian(2, 1) = hessian(1, 2)
     hessian(2, 2) = dot_product(-deriv_c2, -deriv_c2) + dot_product(R, -deriv2_c2)

     ! Bounds Checking
     do i=1, 2
        flag = .False.
        if (pt(i) < low(i)+eps .and. grad(i) >= 0.0) then
           flag = .True.
           pt(i) = low(i)
        end if
        
        if (pt(i) > high(i)-eps .and. grad(i) <= 0.0) then
           flag = .True.
           pt(i) = high(i)
        end if
        
        if (flag) then
           grad(i) = 0.0
           hessian(:, i) = 0.0
           hessian(i, :) = 0.0
           hessian(i, i) = 1.0
        end if
     end do

     ! Zero-cosine check...slightly modifed from The NURBS Book
     if (NORM2(grad) < eps) then
        exit iteration_loop
     end if

     ! Invert the hessian, compute the update and the projected gradient
     call solve_2by2(hessian, grad, update)

     update = -update
     pgrad = dot_product(update, grad)
 
     !Check that this is a descent direction - 
     !otherwise use the negative gradient    
     if ( pgrad >= 0.0 ) then
        update = -grad/NORM2(grad)
        pgrad = dot_product(update, grad)
     end if
     
     step = 1.0
     nDist = 0.0
   
     lineloop: do m=1, nLine
        newpt = pt + step * update
        cflag = .False. ! Check if the constraints are applied
        ! Check if the new point exceeds the bounds
        do i=1, 2
           if (newpt(i) > high(i)) then
              cflag = .True.
              newpt(i) = high(i)
           else if (newpt(i) < low(i)) then
              cflag = .True.
              newpt(i) = low(i)
           end if
        end do
           
        ! Re-evaluate new position
        call eval_curve(newPt(1), t1, k1, coef1, n1, ndim, 1, val1)
        call eval_curve(newPt(2), t2, k2, coef2, n2, ndim, 1, val2)

        ! Check if the new function value is lower, 
        ! otherwise adjust the step size
        R = val1 - val2
        ndist = NORM2(R)
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
     
     if ( m == nLine + 1) then
        dist = ndist
        nLSFail = nLSFail + 1
        if (NLSFail > LSFailMax) then
           exit Iteration_loop
        end if
     else
        NLSFail = 0
        ! Check if there has been no change in the coordinates
        p_diff = NORM2(pt-newpt)
        if (p_diff < eps) then
           exit Iteration_loop
        end if
     end if
     
     pt = newpt
  end do iteration_loop

  s = pt(1)
  t = pt(2)
  diff = R

end subroutine curve_curve

subroutine curve_surface(tc, kc, coefc, tu, tv, ku, kv, coefs, &
     nctlc, nctlu, nctlv, ndim, niter, eps, u, v, s, Diff)

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
  !     coefc   - Real, aray size(ndim, nctlc) coefficients for curve
  !     tu      - Real, Knot vector in u. Length nctlu+ku
  !     tv      - Real, Knot vector in v. Length nctlv+kv
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     coefs   - Real, Array of B-spline coefficients  Size (ndim, nctlv, nctlu)
  !     nctlu   - Integer, Number of control points in u
  !     nctlv   - Integer, Number of control points in v
  !     ndim    - Integer, Spatial Dimension
  !     Niter   - Integer, Maximum number of Netwton iterations
  !     eps     - Real - Parametric convergence measure
  !
  !     Ouput 
  !     u       - Real, u parameter where S(u, v) is closest to C(s)
  !     v       - Real, v parameter where S(u, v) is closest to C(s)
  !     s       - Real, s parameter where S(u, v) is closest to C(s)
  !     diff    - Real Array size(ndim) - Distance between C(s)and S(u, v)

  use precision
  implicit none

  ! Input 
  integer            , intent(in)      :: kc, ku, kv, nctlu, nctlv, nctlc, ndim, niter
  real(kind=realType), intent(in)      :: tc(nctlc+kc), tu(nctlu+ku), tv(nctlv+kv)
  real(kind=realType), intent(in)      :: coefc(ndim, nctlc), coefs(ndim, nctlv, nctlu)
  real(kind=realType), intent(in)      :: eps

  ! Output
  real(kind=realType), intent(inout)   :: u, v, s
  real(kind=realType), intent(out)     :: diff(ndim)

  ! Working
  real(kind=realType)  :: val_s(ndim), deriv_s(ndim, 2), deriv2_s(ndim, 2, 2)
  real(kind=realType)  :: val_c(ndim), deriv_c(ndim), deriv2_C(ndim)
  real(kind=realType)  :: dist, hessian(3,3), grad(3), newpt(3), c
  real(kind=realType)  :: pgrad, update(3), step, p_diff
  real(kind=realType)  :: low(3), high(3), pt(3), R(nDim), ndist, fval, nfval
  integer              :: i, m, ii
  integer              :: NLSFail
  logical              :: flag, cflag

  NLSFail = 0

  ! Set lower and upper bounds of u,v, s based on knot vector
  low(1) = tu(1)
  low(2) = tv(1)
  low(3) = tc(1)
  high(1) = tu(Nctlu+ku)
  high(2) = tv(Nctlv+kv)
  high(3) = tc(Nctlc+kc)

  pt(1) = u
  pt(2) = v
  pt(3) = s

  iteration_loop: do ii=1, niter
     
     ! Evaluate the surface and the curve
     call eval_surface(pt(1), pt(2), tu, tv, ku, kv, coefs, nctlu, nctlv, &
          ndim, 1, 1, val_s)
     call eval_surface_deriv(pt(1), pt(2), tu, tv, ku, kv, coefs, &
          nctlu, nctlv, ndim, deriv_s)
     call eval_surface_deriv2(pt(1), pt(2), tu, tv, ku, kv, coefs, &
          nctlu, nctlv, ndim, deriv2_s)

     call eval_curve(pt(3), tc, kc, coefc, nctlc, ndim, 1, val_c)
     call eval_curve_deriv(pt(3), tc, kc, coefc, nctlc, ndim, deriv_c)
     call eval_curve_deriv2(pt(3), tc, kc, coefc, nctlc, ndim, deriv2_c)

     ! Distance is R, "function value" fval is what we minimize
     R = val_s - val_c
     nDist = NORM2(R)
     fVal = 0.5*nDist**2

     ! Calculate the gradient
     grad(1) = dot_product(R, deriv_s(:,1))
     grad(2) = dot_product(R, deriv_s(:,2))
     grad(3) = dot_product(R, -deriv_c)

     ! Calculate the Hessian
     hessian = 0.0
     hessian(1, 1) = dot_product(deriv_s(:,1), deriv_s(:,1)) + dot_product(R, deriv2_s(:,1,1))
     hessian(1, 2) = dot_product(deriv_s(:,1), deriv_s(:,2)) + dot_product(R, deriv2_s(:,1,2))
     hessian(1, 3) = dot_product(deriv_s(:,1), -deriv_c)

     hessian(2, 1) = hessian(1,2)
     hessian(2, 2) = dot_product(deriv_s(:,2), deriv_s(:,2)) + dot_product(R, deriv2_s(:,2,2))
     hessian(2, 3) = dot_product(deriv_s(:,2), -deriv_c)

     hessian(3, 1) = hessian(1, 3)
     hessian(3, 2) = hessian(2, 3)
     hessian(3, 3) = dot_product(deriv_c, deriv_c) - dot_product(R, deriv2_c)
     
     ! Bounds Checking
     do i=1, 3
        flag = .False.
        if (pt(i) < low(i)+eps .and. grad(i) >= 0.0) then
           flag = .True.
           pt(i) = low(i)
        end if
        
        if (pt(i) > high(i)-eps .and. grad(i) <= 0.0) then
           flag = .True.
           pt(i) = high(i)
        end if
        
        if ( flag ) then
           grad(i) = 0.0
           hessian(:, i) = 0.0
           hessian(i, :) = 0.0
           hessian(i, i) = 1.0
        end if
     end do

     if (NORM2(grad) < eps) then
        exit iteration_loop
     end if

     ! Invert the hessian, compute the update and the projected gradient
     call solve_3by3(hessian, grad, update)

     update = -update
     pgrad = dot_product(update, grad)
 
     !Check that this is a descent direction - 
     !otherwise use the negative gradient    
     if ( pgrad >= 0.0 ) then
        update = -grad/NORM2(grad)
        pgrad = dot_product(update, grad)
     end if
     
     step = 1.0
     nDist = 0.0
   
     lineloop: do m=1, nLine
        newpt = pt + step * update
        cflag = .False. ! Check if the constraints are applied

        ! Check if the new point exceeds the bounds
        do i=1, 3
           if (newpt(i) > high(i)) then
              cflag = .True.
              newpt(i) = high(i)
           else if (newpt(i) < low(i)) then
              cflag = .True.
              newpt(i) = low(i)
           end if
        end do

        ! Evaluate the new point
        call eval_surface(newPt(1), newPt(2), tu, tv, ku, kv, coefs, &
             nctlu, nctlv, ndim, 1, 1, val_s)
        call eval_curve(newPt(3), tc, kc, coefc, nctlc, ndim, 1, val_c)

        ! Distance is R, "function value" fval is what we minimize
        R = val_s - val_c
        nDist = NORM2(R)
        nfVal = 0.5*nDist**2

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
     
     if ( m == nLine+1 ) then
        dist = ndist
        nLSFail = nLSFail + 1

        if (NLSFail > LSFailMax) then ! There is nothing more we can do...
           exit Iteration_loop
        end if
     else
        NLSFail = 0
        ! Check if there has been no change in the coordinates
        p_diff = NORM2(pt-newpt)
        if (p_diff < eps) then
           exit Iteration_loop
        end if
     end if

     pt = newpt
  end do iteration_loop

  u = pt(1)
  v = pt(2)
  s = pt(3)
  diff = R

end subroutine curve_surface

! ------------------------------------------------------------------------------------
!              Globalization Functions
! ------------------------------------------------------------------------------------

subroutine point_curve_start(x0, uu, data, nu, ndim, N, u)
  
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: point_curve_start uses discrete data to determine a good
  !     starting point for the curve projection algorithm

  !     Description of Arguments
  !     Input
  !     x0      - Real, array size (ndim, N), Points we want to globalize
  !     uu      - Real, array size(nu) u-parameter values defining data
  !     data    - Real, array size(ndim, nu) - Data to compare against
  !     nu      - Integer, number of uu data
  !     ndim    - Integer, Spatial dimension
  !     N       - Integer, number of points to check
  !
  !     Ouput
  !     u       - Real, array size(N) - cloested u-parameters

  use precision
  implicit none

  ! Input
  integer              , intent(in) :: nu, ndim, N
  real(kind=realType)  , intent(in) :: x0(ndim, N), uu(nu), data(ndim, nu)
  
  ! Output
  real(kind=realType)  , intent(out) :: u(N)

  ! Working
  real(kind=realType)  :: D
  integer              :: ipt, i

  do ipt=1,N
     D = 1e20
     do i=1,nu
        if (NORM2(X0(:, ipt)-data(:, i)) < D) then
           u(ipt) = uu(i)
           D = NORM2(X0(:, ipt)-data(:, i))
        end if
     end do
  end do

end subroutine point_curve_start

subroutine point_surface_start(x0, uu, vv, data, nu, nv, ndim, N, u, v)
  
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: point_surface_start uses discrete data to determine a good
  !     starting point for the surface projection algorithm

  !     Description of Arguments
  !     Input
  !     x0      - Real, array size (ndim, N), Points we want to globalize
  !     uu      - Real, array size(nu) u-parameter values defining data
  !     vv      - Real, array size(nv) v-parameter values defining data
  !     data    - Real, array size(ndim, nw, nv, nu) - Data to compare against
  !     nu      - Integer, number of uu data
  !     nv      - Integer, number of vv data
  !     ndim    - Integer, Spatial dimension
  !     N       - Integer, number of points to check
  !
  !     Ouput
  !     u       - Real, array size(N) - cloested u-parameters
  !     v       - Real, array size(N) - cloested v-parameters

  use precision
  implicit none

  ! Input
  integer              , intent(in) :: nu, nv, ndim, N
  real(kind=realType)  , intent(in) :: x0(ndim, N), uu(nu), vv(nv)
  real(kind=realType)  , intent(in) :: data(ndim, nv, nu)
  
  ! Output
  real(kind=realType)  , intent(out) :: u(N), v(N)

  ! Working
  real(kind=realType)  :: D
  integer              :: ipt, i, j
  
  do ipt=1,N
     D = 1e20
     do i=1,nu
        do j=1,nv
           if (NORM2(X0(:, ipt)-data(:, j, i)) < D) then
              u(ipt) = uu(i)
              v(ipt) = vv(j)
              D = NORM2(X0(:, ipt)-data(:, j, i))
           end if
        end do
     end do
  end do

end subroutine point_surface_start

subroutine point_volume_start(x0, uu, vv, ww, data, nu, nv, nw, ndim, N, u, v, w)
  
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: point_volume_start uses discrete data to determine a good
  !     starting point for the volume projection algorithm

  !     Description of Arguments
  !     Input
  !     x0      - Real, array size (ndim, N), Points we want to globalize
  !     uu      - Real, array size(nu) u-parameter values defining data
  !     vv      - Real, array size(nv) v-parameter values defining data
  !     ww      - Real, array size(nw) w-parameter values defining data
  !     data    - Real, array size(ndim, nw, nv, nu) - Data to compare against
  !     nu      - Integer, number of uu data
  !     nv      - Integer, number of vv data
  !     nw      - Integer, number of ww data
  !     ndim    - Integer, Spatial dimension
  !     N       - Integer, number of points to check
  !
  !     Ouput
  !     u       - Real, array size(N) - cloested u-parameters
  !     v       - Real, array size(N) - cloested v-parameters
  !     w       - Real, array size(N) - cloested w-parameters

  use precision
  implicit none

  ! Input
  integer              , intent(in) :: nu, nv, nw, ndim, N
  real(kind=realType)  , intent(in) :: x0(ndim, N), uu(nu), vv(nv), ww(nw)
  real(kind=realType)  , intent(in) :: data(ndim, nw, nv, nu)
  
  ! Output
  real(kind=realType)  , intent(out) :: u(N), v(N), w(N)

  ! Working
  real(kind=realType)  :: D
  integer              :: ipt, i, j, k

  do ipt=1,N
     D = 1e20
     do i=1,nu
        do j=1,nv
           do k=1,nw
              if (NORM2(X0(:, ipt)-data(:,k, j, i)) < D) then
                 u(ipt) = uu(i)
                 v(ipt) = vv(j)
                 w(ipt) = ww(k)
                 D = NORM2(X0(:, ipt)-data(:,k, j, i))
              end if
           end do
        end do
     end do
  end do

end subroutine point_volume_start

subroutine curve_curve_start(data1, uu1, data2, uu2, nu1, nu2, ndim, s1, s2)
  
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: point_surface_start uses discrete data to determine a good
  !     starting point for the surface projection algorithm

  !     Description of Arguments
  !     Input
  !     data1   - Real, array size (ndim, nu1), Points from first curve
  !     uu1     - Real, array size (nu1), parameter values from first curve
  !     data2   - Real, array size (ndim, nu2), Points from second curve
  !     uu2     - Real, array size (nu2), parameter values from second curve
  !     nu1     - Integer, number of points in data1
  !     nu2     - Integer, number of points in data2
  !     ndim    - Integer, Spatial dimension
  !
  !     Ouput
  !     s1      - Real, parameter value on curve1
  !     s2      - Real, parameter value on curve2

  use precision
  implicit none

  ! Input
  integer              , intent(in) :: nu1, nu2, ndim
  real(kind=realType)  , intent(in) :: data1(ndim, nu1), uu1(nu1), data2(ndim, nu2), uu2(nu2)
  
  ! Output
  real(kind=realType)  , intent(out) :: s1, s2

  ! Working
  real(kind=realType)  :: D
  integer              :: i, j

  D = 1e20
  do i=1,nu1
     do j=1,nu2
        if (NORM2(data1(:,i) - data2(:,j)) < D) then
           s1 = uu1(i)
           s2 = uu2(j)
           D = NORM2(data1(:,i) - data2(:,j))
        end if
     end do
  end do

end subroutine curve_curve_start

subroutine curve_surface_start(data1, uu1, data2, uu2, vv2, nu1, nu2, nv2, ndim, s, u, v)
  
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: point_surface_start uses discrete data to determine a good
  !     starting point for the surface projection algorithm

  !     Description of Arguments
  !     Input
  !     data1   - Real, array size (ndim, nu1), Points from curve
  !     uu1     - Real, array size (nu1), parameter values from curve
  !     data2   - Real, array size (ndim, nv2, nu2), Points from surface
  !     uu2     - Real, array size (nu2), u parameter values from surface
  !     vv2     - Real, array size (nv2), v parameter values from surface
  !     nu1     - Integer, number of points in data1
  !     nu2     - Integer, number of u-points in data2
  !     nv2     - Integer, number of v-points in data2
  !     ndim    - Integer, Spatial dimension
  !
  !     Ouput
  !     s       - Real, parameter value on curve
  !     u       - Real, u-parameter value on surface
  !     v       - Real, v-parameter value on surface

  use precision
  implicit none

  ! Input
  integer              , intent(in) :: nu1, nu2, nv2, ndim
  real(kind=realType)  , intent(in) :: data1(ndim, nu1), uu1(nu1)
  real(kind=realtype)  , intent(in) :: data2(ndim, nv2, nu2), uu2(nu2), vv2(nv2)
  
  ! Output
  real(kind=realType)  , intent(out) :: s, u, v

  ! Working
  real(kind=realType)  :: D
  integer              :: i, j, k

  D = 1e20
  do i=1,nu1
     do j=1,nu2
        do k=1,nv2
           if (NORM2(data1(:,i) - data2(:, k, j)) < D) then
              s = uu1(i)
              u = uu2(j)
              v = vv2(k)
              D = NORM2(data1(:,i) - data2(:, k, j))
           end if
        end do
     end do
  end do
end subroutine curve_surface_start

subroutine solve_2by2(A, b, x)
  
  use precision
  implicit none

  ! Solve a 2 x 2 system  -- With NO checking
  real(kind=realType), intent(in) :: A(2, 2), b(2)
  real(kind=realType), intent(out) :: x(2)
  real(kind=realType)         :: idet, det

  det = A(1, 1)*A(2, 2)-A(1, 2)*A(2, 1)
  if (det == 0) then
     X = B
  else
     idet = 1.0/det
     X(1) = idet*(B(1)*A(2, 2) - B(2)*A(1, 2))
     X(2) = idet*(A(1, 1)*B(2) - B(1)*A(2, 1))
  end if

end subroutine solve_2by2

subroutine solve_3by3(A, b, x)

  use precision
  implicit none

  ! Solve a 3 x 3 system  -- With NO checking
  real(kind=realType), intent(in) :: A(3, 3), b(3)
  real(kind=realType), intent(out) :: x(3)
  real(kind=realType)         :: idet

  idet = 1/(A(1,1)*(A(3,3)*A(2,2)-A(3,2)*A(2,3))-A(2,1)*(A(3,3)*A(1,2)-A(3,2)*A(1,3))+A(3,1)*(A(2,3)*A(1,2)-A(2,2)*A(1,3)))
  x(1) = idet*( b(1)*(A(3,3)*A(2,2)-A(3,2)*A(2,3)) - b(2)*(A(3,3)*A(1,2)-A(3,2)*A(1,3)) + b(3)*(A(2,3)*A(1,2)-A(2,2)*A(1,3)))
  x(2) = idet*(-b(1)*(A(3,3)*A(2,1)-A(3,1)*A(2,3)) + b(2)*(A(3,3)*A(1,1)-A(3,1)*A(1,3)) - b(3)*(A(2,3)*A(1,1)-A(2,1)*A(1,3)))
  x(3) = idet*( b(1)*(A(3,2)*A(2,1)-A(3,1)*A(2,2)) - b(2)*(A(3,2)*A(1,1)-A(3,1)*A(1,2)) + b(3)*(A(2,2)*A(1,1)-A(2,1)*A(1,2)))


  ! | a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
  ! | a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
  ! | a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |

  ! DET  =  a11(a33a22-a32a23)-a21(a33a12-a32a13)+a31(a23a12-a22a13)

end subroutine solve_3by3

subroutine line_plane(ia, vc, p0, v1, v2, n, sol, pid, n_sol)

  ! Check a line against multiple planes
  !
  ! ia:   The initial point
  ! vc:   The search vector from the initial point
  ! p0:   Vectors to the triangle origins
  ! v1:   Vector along the first triangle direction
  ! v2:   Vector along the second triangle direction
  !  n:   Number of triangles to search
  ! sol:  Solution vector - parametric positions + physical coordiantes
  ! nsol: Number of solutions
  !
  ! Solve for the scalars: alpha, beta, gamma such that:
  !    ia + alpha*vc = p0 + beta*v1 + gamma*v2
  !    ia - p0 = [ - vc ; v1 ; v2 ][ alpha ]
  !                                [ beta  ]
  !                                [ gamma ]
  !
  ! alpha >= 0: The point lies above the initial point
  ! alpha  < 0: The point lies below the initial point
  !
  ! The domain of the triangle is defined by: 
  !     beta + gamma = 1
  ! and 
  !     0 < beta, gamma < 1

  use precision
  implicit none 

  ! Input
  integer, intent(in) :: n
  real(kind=realType) , intent(in) :: ia(3), vc(3), p0(3, n), v1(3, n), v2(3, n)

  ! Output
  integer, intent(out) :: n_sol
  real(kind=realType), intent(out) :: sol(6, n)
  integer, intent(out) :: pid(n)

  ! Worling 
  integer :: i
  real(kind=realType) :: A(3, 3), rhs(3), x(3)
  
  A(:, 1) = -vc
  n_sol = 0
  sol(:, :) = 0.0
  do i=1, n
     A(:, 2) = v1(:, i)
     A(:, 3) = v2(:, i)
     rhs = ia-p0(:, i)
     
     call solve_3by3(A, rhs, x)
     
     ! if (x(1) .ge. 0.00 .and. x(1) .le. 1.00 .and. &
     !     x(2) .ge. 0.00 .and. x(2) .le. 1.00 .and. &
     !     x(3) .ge. 0.00 .and. x(3) .le. 1.00 .and. &
     !     x(2)+x(3) .le. 1.00) then

     if  (x(2) .ge. 0.0 .and. x(2) .le. 1.0 .and. &
          x(3) .ge. 0.0 .and. x(3) .le. 1.0 .and. &
          x(2) + x(3) .le. 1.0) then

        n_sol = n_sol + 1
        sol(1:3, n_sol) = x  ! t, u, v parametric locations
        sol(4:6, n_sol) = ia + x(1)*vc ! Actual point value
        pid(n_sol) = i
     end if
  end do
end subroutine line_plane

subroutine plane_line(ia, vc, p0, v1, v2, n, sol, n_sol)

  ! Check a plane against multiple lines

  use precision 
  implicit none 

  ! Input
  integer, intent(in) :: n
  real(kind=realType) , intent(in) :: ia(3, n), vc(3, n), p0(3), v1(3), v2(3)

  ! Output
  integer, intent(out) :: n_sol
  real(kind=realType), intent(out) :: sol(6, n)

  ! Worling 
  integer :: i, ind(n)
  real(kind=realType) :: A(3, 3), rhs(3), x(3)

  n_sol = 0
  sol(:, :) = 0.0
  A(:, 2) = v1(:)
  A(:, 3) = v2(:)
  
  do i=1, n
  
     A(:, 1) = -vc(:, i)
     rhs = ia(:, i)-p0(:)
     
     call solve_3by3(A, rhs, x)
     
     if (x(1) .ge. 0.00 .and. x(1) .le. 1.00 .and. &
         x(2) .ge. 0.00 .and. x(2) .le. 1.00 .and. &
         x(3) .ge. 0.00 .and. x(3) .le. 1.00 .and. &
         x(2)+x(3) .le. 1.00) then

        n_sol = n_sol + 1
        sol(1:3, n_sol) = x  ! t, u, v parametric locations
        sol(4:6, n_sol) = ia(:, i) + x(1)*vc(:, i) ! Actual point value
        ind(n_sol) = i

     end if
  end do
end subroutine plane_line

subroutine point_plane(pt, p0, v1, v2, n, sol, n_sol, best_sol)

  use precision
  implicit none
 
  ! Input
  integer, intent(in) :: n
  real(kind=realType) , intent(in) :: pt(3), p0(3, n), v1(3, n), v2(3, n)

  ! Output
  integer, intent(out) :: n_sol, best_sol
  real(kind=realType), intent(out) :: sol(6, n)

  ! Working 
  integer :: i, ind(n)
  real(kind=realType) :: A(2, 2), rhs(2), x(2), r(3), D, D0
  

  n_sol = 0
  sol(:, :) = 0.0
  do i=1, n
     A(1, 1) = v1(1, i)**2 + v1(2, i)**2 + v1(3, i)**2
     A(1, 2) = v1(1, i)*v2(1, i) + v1(2, i)*v2(2, i) + v1(3, i)+v2(3, i)
     A(2, 1) = A(1, 2)
     A(2, 2) = v2(1, i)**2 + v2(2, i)**2 + v2(3, i)**2
     r = p0(:, i)-pt
     rhs(1) = r(1)*v1(1, i) + r(2)*v1(2, i) + r(3)*v1(3, i)
     rhs(2) = r(1)*v2(1, i) + r(2)*v2(2, i) + r(3)*v2(3, i)
     
     call solve_2by2(A, rhs, x)
     
     if (x(1) .ge. 0.00 .and. x(1) .le. 1.00 .and. &
         x(2) .ge. 0.00 .and. x(2) .le. 1.00 .and. &
         x(1)+x(2) .le. 1.00) then

        n_sol = n_sol + 1
        sol(2:3, n_sol) = x  ! t, u, v parametric locations
        sol(4:6, n_sol) = 0.0 ! Actual point value
        ind(n_sol) = i
     end if
  end do

  ! Now post-process to get the closest one
  best_sol = 1
  D0 = NORM2(p0(:, ind(1)) + sol(2, ind(1))*v1(:, ind(1)) + sol(3, ind(1))*v2(:, ind(1)))

  do i=1, n_sol
     D = NORM2(p0(:, ind(i)) + sol(2, ind(i))*v1(:, ind(i)) + sol(3, ind(i))*v2(:, ind(i)))
     if (D<D0) then
        D0 = D
        best_sol = i
     end if
  end do
end subroutine point_plane

function dotproduct(x1, x2, n)
  use precision 
  implicit none

  integer, intent(in) :: n
  integer :: i
  real(kind=realType), intent(in) :: x1(n), x2(n)

  real(kind=realType) ::  dotproduct 
  dotproduct = 0.0
  do i=1, n
     dotproduct = dotproduct + x1(i)*x2(i)
  end do

end function dotproduct

