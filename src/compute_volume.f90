subroutine volume_jacobian_wrap(u, v, w, tu, tv, tw, ku, kv, kw, nctlu, nctlv, nctlw, &
     nu, nv, nw, vals, row_ptr, col_ind)
  use precision
  implicit none

  ! Input
  integer            , intent(in)      :: ku, kv, kw, nctlu, nctlv, nctlw, nu, nv, nw
  real(kind=realType), intent(in)      :: u(nw, nv, nu), v(nw, nv, nu), w(nw, nv, nu)
  real(kind=realType), intent(in)      :: tu(nctlu+ku), tv(nctlv+kv), tw(nctlw+kw)

  ! Output
  real(kind=realType), intent(out)     :: vals(nu*nv*nw*ku*kv*kw)
  integer            , intent(out)     :: col_ind(nu*nv*nw*ku*kv*kw), row_ptr(nu*nv*nw+1)

  ! Working
  real(kind=realType)                  :: basisu(ku), basisv(kv), basisw(kw)
  integer                              :: i, j, k, ii, jj, kk, counter, c1, c2
  integer                              :: ileftu, ileftv, ileftw

  counter = 1
  do i=1, nu
     do j = 1, nv
        do k=1, nw
           ! U 
           call findSpan(u(k, j, i), ku, tu, nctlu, ileftu)
           call basis(tu, nctlu, ku, u(k, j, i), ileftu, basisu)

           ! V
           call findSpan(v(k, j, i), kv, tv, nctlv, ileftv)
           call basis(tv, nctlv, kv, v(k, j, i), ileftv, basisv)

           ! W
           call findSpan(w(k, j, i), kw, tw, nctlw, ileftw)
           call basis(tw, nctlw, kw, w(k, j, i), ileftw, basisw)

           row_ptr( (i-1)*nv*nw + (j-1)*nw + k  ) = counter-1
           do ii=1, ku
              c1 = (ileftu-ku+ii-1)*Nctlv*Nctlw
              do jj = 1, kv
                 c2 = (ileftv-kv+jj-1)*Nctlw
                 do kk=1, kw
                    col_ind(counter) = c1 + c2 + (ileftw-kw+kk-1)
                    vals(counter) = basisu(ii)*basisv(jj)*basisw(kk)
                    counter = counter + 1
                 end do
              end do
           end do
        end do
     end do
  end do
  row_ptr(nu*nv*nw+1) = counter-1

end subroutine volume_jacobian_wrap
