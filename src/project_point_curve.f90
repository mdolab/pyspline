subroutine project_point_curve(x0,t,k,coef,nctl,ndim,Niter,eps1,eps2,s,Diff)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: project_point attempts to solve the point inversion problem
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
  double precision                      :: val0(ndim),s0,fuck_val
  integer                               :: i,j
  double precision                      :: D,D0
  integer                               :: n ! Huristic Value

  ! Functions
  double precision                      :: norm

  ! Is 3  Good Here?? Huristic!
  n=3
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
     
     if ( norm(dot_product(deriv,Diff),ndim)/ (norm(deriv,ndim)*norm(Diff,ndim)) <= eps2) then
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

end subroutine project_point_curve


