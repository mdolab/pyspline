subroutine surface_jacobian_wrap(u, v, tu, tv, ku, kv, nctlu, nctlv, nu, nv, vals, &
     row_ptr, col_ind)

  use precision
  implicit none

  ! Input
  integer            , intent(in)      :: ku, kv, nctlu, nctlv, nu, nv
  real(kind=realType), intent(in)      :: u(nv, nu), v(nv, nu)
  real(kind=realType), intent(in)      :: tu(nctlu+ku), tv(nctlv+kv)

  ! Output
  real(kind=realType), intent(out)     :: vals(nu*nv*ku*kv)
  integer            , intent(out)     :: col_ind(nu*nv*ku*kv), row_ptr(nu*nv+1)

  ! Working
  real(kind=realType)                  :: basisu(ku), basisv(kv)
  integer                              :: ileftu, ileftv
  integer                              :: i, j, ii, jj, counter

  counter = 1
  do i=1, nu
     do j = 1, nv
        ! Get u interval
        call findSpan(u(j, i), ku, tu, nctlu, ileftu)
        call basis(tu, nctlu, ku, u(j, i), ileftu, basisu)
        
        ! Get v interval
        call findSpan(v(j, i), kv, tv, nctlv, ileftv)
        call basis(tv, nctlv, kv, v(j, i), ileftv, basisv)

        row_ptr( (i-1)*nv + j  ) = counter-1
        do ii=1, ku
           do jj = 1, kv
              col_ind(counter) = (ileftu-ku+ii-1)*Nctlv + (ileftv-kv+jj-1)
              vals(counter) = basisu(ii)*basisv(jj)
              counter = counter + 1
              
           end do
        end do
     end do
  end do
  row_ptr(nu*nv+1) = counter-1
  
end subroutine surface_jacobian_wrap

subroutine surface_para_corr(tu, tv, ku, kv, u, v, coef, nctlu, nctlv, ndim, nu, nv, X, rms)

  ! Do Hoschek parameter correction
  use precision
  implicit none

  ! Input/Output
  integer              , intent(in)      :: ku, kv, nctlu, nctlv, ndim, nu, nv
  real(kind=realType)  , intent(in)      :: tu(ku+nctlu), tv(kv+nctlv)
  real(kind=realType)  , intent(inout)   :: u(nv, nu), v(nv, nu)
  real(kind=realType)  , intent(in)      :: coef(ndim, nctlv, nctlu)
  real(kind=realType)  , intent(in)      :: X(ndim, nv, nu)
  real(kind=realType)  , intent(out)     :: rms

  ! Working
  integer                               :: i, j, jj, max_inner_iter
  real(kind=realType)                   :: D(ndim), D2(ndim)
  real(kind=realType)                   :: val(ndim), deriv(ndim, 2), deriv2(ndim, 2, 2)
  real(kind=realType)                   :: u_tilde, v_tilde
  real(kind=realType)                   :: A(2, 2), ki(2), delta(2)

  max_inner_iter = 10
  rms = 0.0

  do i=2, nu-1
     do j = 2, nv-2
        call eval_surface(u(j, i), v(j, i), tu, tv, ku, kv, coef, nctlu, nctlv, ndim, val)
        call eval_surface_deriv(u(j, i), v(j, i), tu, tv, ku, kv, coef, nctlu, nctlv, ndim, deriv)
        call eval_surface_deriv2(u(j, i), v(j, i), tu, tv, ku, kv, coef, nctlu, nctlv, ndim, deriv2)

        D = val-X(:, j, i)

        A(1, 1) = NORM2(deriv(:, i))**2 + dot_product(D, deriv2(:, 1, 1))
        A(1, 2) = dot_product(deriv(:, 1), deriv(:, 2)) + dot_product(D, deriv2(:, 1, 2))
        A(2, 1) = A(1, 2)
        A(2, 2) = NORM2(deriv(:, 2))**2 + dot_product(D, deriv2(:, 2, 2))
        
        ki(1) = -dot_product(D, deriv(:, 1))
        ki(2) = -dot_product(D, deriv(:, 2))
        
        call solve_2by2(A, ki, delta)
        
        if (j .eq.1 .or. j .eq. nv) then
           delta(1) = 0.0
        end if
        if (i .eq.1 .or. i .eq. nu) then
           delta(2) = 0.0
        end if
        inner_loop: do jj=1, max_inner_iter
           u_tilde = u(j, i) + delta(1)
           v_tilde = v(j, i) + delta(2)

           call eval_surface(u_tilde, v_tilde, tu, tv, ku, kv, coef, nctlu, nctlv, ndim, val)
           D2 = val-X(i, j, :)
           if (NORM2(D) .ge. NORM2(D2)) then
              u(j, i) = u_tilde
              v(j, i) = v_tilde
              exit inner_loop
           else
              delta = delta*0.5
           end if
        end do inner_loop
     end do
  end do

  ! Lets redo the full RMS
  rms = 0.0

  do i=1, nu
     do j=1, nv
        call eval_surface(u(j, i), v(j, i), tu, tv, ku, kv, coef, nctlu, nctlv, ndim, val)
        D = X(i, j, :)-val
        rms = rms + dot_product(D, D)
     end do
  end do
  rms = sqrt(rms/(nu*nv))

end subroutine surface_para_corr

function compute_rms_surface(tu, tv, ku, kv, u, v, coef, nctlu, nctlv, ndim, nu, nv, X)
 ! Do Hoschek parameter correction
  use precision
  implicit none

  ! Input/Output
  integer              , intent(in)      :: ku, kv, nctlu, nctlv, ndim, nu, nv
  real(kind=realType)  , intent(in)      :: tu(ku+nctlu), tv(kv+nctlv)
  real(kind=realType)  , intent(inout)   :: u(nv, nu), v(nv, nu)
  real(kind=realType)  , intent(in)      :: coef(ndim, nctlv, nctlu)
  real(kind=realType)  , intent(in)      :: X(ndim, nv, nu)
 
  ! Working
  integer                            :: i, j, idim
  real(kind=realType)                   :: val(ndim), D(ndim)
  real(kind=realType)                   :: compute_rms_surface

  compute_rms_surface = 0.0
  do i=1, nu
     do j=1, nv
        call eval_surface(u(j, i), v(j, i), tu, tv, ku, kv, coef, nctlu, nctlv, ndim, val)
        D = val-X(:, j, i)
        do idim=1, ndim
           compute_rms_surface = compute_rms_surface + D(idim)**2
        end do
     end do
  end do
  compute_rms_surface = sqrt(compute_rms_surface/(nu*nv))

end function compute_rms_surface

