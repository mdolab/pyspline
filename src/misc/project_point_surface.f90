subroutine project_point_surface(x0,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,niter,eps1,eps2,u,v,Diff)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: project_point_surface attempts to solve the point inversion problem for a surface
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
  double precision, intent(out)         :: u,v,diff(ndim)


  ! Working
  double precision                      :: val(ndim),deriv(2,ndim),deriv2(2,2,ndim)
  double precision                      :: val0(ndim),s0,fuck_val
  integer                               :: i,j,ii,jj
  double precision                      :: D,D0,u0,v0,delta(2)
  double precision                      :: A(2,2),ki(2)
  integer                               :: n ! Huristic Value

  ! Functions     
  double precision                      :: norm

  n = 3

  ! First we will evaluate the surface at n points inside each knot span in each direction

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

  ! Now we have u0 and s0 so we can do the newton 

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
end subroutine project_point_surface


subroutine solve_2by2(A,b,x)
  ! Solve a 2 x 2 system  -- With NO checking
  double precision         :: A(2,2),b(2),x(2)
  double precision         :: idet

  idet = 1/(A(1,1)*A(2,2)-A(1,2)*A(2,1))
  
  x(1) = idet*(A(2,2)*b(1) - A(1,2)*b(2))
  x(2) = idet*(-A(2,1)*b(1) + A(1,1)*b(2))

end subroutine solve_2by2

  
