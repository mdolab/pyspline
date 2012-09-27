function bvalu(t, coef, nctl, K, ideriv, s)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract bvalu is the generic 1D spline interpolant
  !              function. It is similar to the function of the same
  !              name from the tensbs library without the work arrays
  !
  !     Description of Arguments
  !     Input
  !     t       - Real, Knot vector. Length nctl+k
  !     coef    - Real, Array of b-sline coefficients  Size (nctl)
  !     nctl    - Integer, Number of control points
  !     k       - Integer, order of B-spline 
  !     ideriv  - Integer, degree of derivatve to evaluate
  !     s       - Real, Position to evaluate 
  !
  !     Ouput bvalu - Real, Evaluated point or derivative

  use precision
  implicit none

  ! Input
  integer         , intent(in)     :: k, nctl, ideriv
  real(kind=realType), intent(in)  :: s
  real(kind=realType), intent(in)  :: t(nctl+k)
  real(kind=realType), intent(in)  :: coef(nctl)

  ! Output
  real(kind=realType)              :: bvalu

  ! Working
  integer                          :: i, l, idim, istart, ileft
  real(kind=realType)              :: B(k), Bd(k,k)

  bvalu = 0.0

  ! Find knot span
  call findSpan(s, k, t, nctl, ileft)
  istart = ileft-k
  if (ideriv ==0) then ! Use basis function
     call basis(t, nctl, k, s, ileft, B)
     do l=1, k
        bvalu = bvalu + B(l)*coef(istart+l)
     end do
  else
     ! Evalute derivative basis functions up to depree specified by
     ! ideriv
     call derivBasis(t, nctl, k, s, ileft, ideriv, Bd)
     do l=1, k
        bvalu = bvalu + Bd(ideriv+1,l)*coef(istart+l)
     end do
  end if
end function bvalu

function b2val(u, v, idu, idv, tu, tv, nctlu, nctlv, ku, kv, coef)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract bvalu is the generic 2D spline interpolant
  !              function. It is similar to the function of the same
  !              name from the tensbs library without the work array
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, parametric position 1
  !     v       - Real, parametric position 2
  !     idu     - Integer, derivative in u to evaluate
  !     idv     - Integer, derivative in v to evaluate
  !     tu      - Real, U Knot vector. Length nctlu+ku
  !     tv      - Real, V Knot vector. Length nctlv+kv
  !     Nctlu   - Integer, Number of u control points
  !     Nctlv   - Integer, Number of v control points
  !     ku      - Integer, order of spline in u direction
  !     kv      - Integer, order of spline in v direction
  !     coef    - Real, Array size Nctlu x Nctlv, B-spline coefficients
  !
  !     Ouput 
  !     b2val   - Real, Evaluated point or derivative

  use precision
  implicit none

  ! Input
  integer         , intent(in)     :: ku, kv, nctlu, nctlv, idu, idv
  real(kind=realType), intent(in)  :: u,v, tu(nctlu+ku), tv(nctlv+kv)
  real(kind=realType), intent(in)  :: coef(nctlu,nctlv)

  ! Output
  real(kind=realType)              :: b2val

  ! Working
  integer                          :: i, j
  integer                          :: ileftu, ileftv, istartu, istartv
  real(kind=realType)              :: Bu(ku), Bv(kv)
  real(kind=realType)              :: Bud(ku, ku), Bvd(kv, kv)

  ! Find knot spans
  call findSpan(u, ku, tu, nctlu, ileftu)
  istartu = ileftu-ku

  call findSpan(v, kv, tv, nctlv, ileftv)
  istartv = ileftv-kv

  b2val = 0.0

  ! Check if we can just use regular basis functions:
  if (idu + idv == 0) then
     call basis(tu, nctlu, ku, u, ileftu, Bu)
     call basis(tv, nctlv, kv, v, ileftv, Bv)

     do i=1, ku
        do j=1, kv
           b2val = b2val + Bu(i)*Bv(j)*coef(istartu+i, istartv+j)
        end do
     end do
  else
     ! We need derivative basis functions
     call derivBasis(tu, nctlu, ku, u, ileftu, idu, Bu)
     call derivBasis(tv, nctlv, kv, v, ileftv, idv, BV)

     do i=1, ku
        do j=1, kv
           b2val = b2val + Bud(idu+1,i)*Bvd(idv+1,j)*coef(istartu+i, istartv+j)
        end do
     end do
  end if

end function b2val

function b3val(u, v, w, idu, idv, idw, tu, tv, tw, nctlu, nctlv, nctlw, &
     ku, kv, kw, coef)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract bvalu is the generic 2D spline interpolant
  !              function. It is similar to the function of the same
  !              name from the tensbs library without the work array
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, parametric position 1
  !     v       - Real, parametric position 2
  !     w       - Real, parametric position 3
  !     idu     - Integer, derivative in u to evaluate
  !     idv     - Integer, derivative in v to evaluate
  !     idw     - Integer, derivative in w to evaluate
  !     tu      - Real, U Knot vector. Length nctlu+ku
  !     tv      - Real, V Knot vector. Length nctlv+kv
  !     tv      - Real, W Knot vector. Length nctlw+kw
  !     Nctlu   - Integer, Number of u control points
  !     Nctlv   - Integer, Number of v control points
  !     Nctlw   - Integer, Number of w control points
  !     ku      - Integer, order of spline in u direction
  !     kv      - Integer, order of spline in v direction
  !     kw      - Integer, order of spline in w direction
  !     coef    - Real, Array size Nctlu x Nctlv x Nctlw, B-spline coefficients
  !
  !     Ouput 
  !     b3val   - Real, Evaluated point or derivative

  use precision
  implicit none

  ! Input
  integer         , intent(in)     :: ku, kv, kw, nctlu, nctlv, nctlw
  integer         , intent(in)     :: idu, idv, idw
  real(kind=realType), intent(in)  :: u, v, w
  real(kind=realType), intent(in)  :: tu(nctlu+ku), tv(nctlv+kv), tw(nctlw+kw)
  real(kind=realType), intent(in)  :: coef(nctlu,nctlv,nctlw)

  ! Output
  real(kind=realType)              :: b3val

  ! Working
  integer                          :: i, j, k
  integer                          :: ileftu, ileftv, ileftw
  integer                          :: istartu, istartv, istartw
  real(kind=realType)              :: Bu(ku), Bv(kv), Bw(kw)
  real(kind=realType)              :: Bud(ku, ku), Bvd(kv, kv), Bwd(kw, kw)

  ! Find knot spans
  call findSpan(u, ku, tu, nctlu, ileftu)
  istartu = ileftu-ku

  call findSpan(v, kv, tv, nctlv, ileftv)
  istartv = ileftv-kv

  call findSpan(w, kw, tw, nctlw, ileftw)
  istartw = ileftw-kw

  b3val = 0.0

  ! Check if we can just use regular basis functions:
  if (idu + idv + idw == 0) then
     call basis(tu, nctlu, ku, u, ileftu, Bu)
     call basis(tv, nctlv, kv, v, ileftv, Bv)
     call basis(tw, nctlw, kw, w, ileftw, Bw)

     do i=1, ku
        do j=1, kv
           do k=1, kw
              b3val = b3val + Bu(i)*Bv(j)*Bw(k)*&
                   coef(istartu+i, istartv+j, istartw+k)
           end do
        end do
     end do
  else
     ! We need derivative basis functions
     call derivBasis(tu, nctlu, ku, u, ileftu, idu, Bu)
     call derivBasis(tv, nctlv, kv, v, ileftv, idv, Bv)
     call derivBasis(tw, nctlw, kw, w, ileftw, idw, Bw)

     do i=1, ku
        do j=1, kv
           do k=1, kw
              b3val = b3val + Bud(idu+1, i)*Bvd(idv+1, j)*Bwd(idw+1, k) * &
                   coef(istartu+i, istartv+j, istartw+k)
           end do
        end do
     end do
  end if

end function b3val

! Also define complex versions --- these are only here for
! TACS. pySpline doesn't need them

function cbvalu(t, coef, nctl, K, ideriv, s)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract bvalu is the generic 1D spline interpolant
  !              function. It is similar to the function of the same
  !              name from the tensbs library without the work arrays
  !
  !     Description of Arguments
  !     Input
  !     t       - Real, Knot vector. Length nctl+k
  !     coef    - Complex, Array of b-sline coefficients  Size (nctl)
  !     nctl    - Integer, Number of control points
  !     k       - Integer, order of B-spline 
  !     ideriv  - Integer, degree of derivatve to evaluate
  !     s       - Real, Position to evaluate 
  !
  !     Ouput bvalu - Complex, Evaluated point or derivative

  use precision
  implicit none

  ! Input
  integer         , intent(in)        :: k, nctl, ideriv
  real(kind=realType), intent(in)     :: s
  real(kind=realType), intent(in)     :: t(nctl+k)
  complex(kind=realType), intent(in)  :: coef(nctl)

  ! Output
  complex(kind=realType)              :: cbvalu

  ! Working
  integer                          :: i, l, idim, istart, ileft
  real(kind=realType)              :: B(k), Bd(k,k)

  cbvalu = cmplx(0.0, 0.0)

  ! Find knot span
  call findSpan(s, k, t, nctl, ileft)
  istart = ileft-k
  if (ideriv ==0) then ! Use basis function
     call basis(t, nctl, k, s, ileft, B)
     do l=1, k
        cbvalu = cbvalu + B(l)*coef(istart+l)
     end do
  else
     ! Evalute derivative basis functions up to depree specified by
     ! ideriv
     call derivBasis(t, nctl, k, s, ileft, ideriv, Bd)
     do l=1, k
        cbvalu = cbvalu + Bd(ideriv+1,l)*coef(istart+l)
     end do
  end if
end function cbvalu

function cb2val(u, v, idu, idv, tu, tv, nctlu, nctlv, ku, kv, coef)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract bvalu is the generic 2D spline interpolant
  !              function. It is similar to the function of the same
  !              name from the tensbs library without the work array
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, parametric position 1
  !     v       - Real, parametric position 2
  !     idu     - Integer, derivative in u to evaluate
  !     idv     - Integer, derivative in v to evaluate
  !     tu      - Real, U Knot vector. Length nctlu+ku
  !     tv      - Real, V Knot vector. Length nctlv+kv
  !     Nctlu   - Integer, Number of u control points
  !     Nctlv   - Integer, Number of v control points
  !     ku      - Integer, order of spline in u direction
  !     kv      - Integer, order of spline in v direction
  !     coef    - Complex, Array size Nctlu x Nctlv, B-spline coefficients
  !
  !     Ouput 
  !     b2val   - Complex, Evaluated point or derivative

  use precision
  implicit none

  ! Input
  integer         , intent(in)        :: ku, kv, nctlu, nctlv, idu, idv
  real(kind=realType), intent(in)     :: u,v, tu(nctlu+ku), tv(nctlv+kv)
  complex(kind=realType), intent(in)  :: coef(nctlu,nctlv)

  ! Output
  complex(kind=realType)              :: cb2val

  ! Working
  integer                          :: i, j
  integer                          :: ileftu, ileftv, istartu, istartv
  real(kind=realType)              :: Bu(ku), Bv(kv)
  real(kind=realType)              :: Bud(ku, ku), Bvd(kv, kv)

  ! Find knot spans
  call findSpan(u, ku, tu, nctlu, ileftu)
  istartu = ileftu-ku

  call findSpan(v, kv, tv, nctlv, ileftv)
  istartv = ileftv-kv

  cb2val = cmplx(0.0, 0.0)

  ! Check if we can just use regular basis functions:
  if (idu + idv == 0) then
     call basis(tu, nctlu, ku, u, ileftu, Bu)
     call basis(tv, nctlv, kv, v, ileftv, Bv)

     do i=1, ku
        do j=1, kv
           cb2val = cb2val + Bu(i)*Bv(j)*coef(istartu+i, istartv+j)
        end do
     end do
  else
     ! We need derivative basis functions
     call derivBasis(tu, nctlu, ku, u, ileftu, idu, Bu)
     call derivBasis(tv, nctlv, kv, v, ileftv, idv, BV)

     do i=1, ku
        do j=1, kv
           cb2val = cb2val + Bud(idu+1,i)*Bvd(idv+1,j)*coef(istartu+i, istartv+j)
        end do
     end do
  end if

end function cb2val

function cb3val(u, v, w, idu, idv, idw, tu, tv, tw, nctlu, nctlv, nctlw, &
     ku, kv, kw, coef)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract bvalu is the generic 2D spline interpolant
  !              function. It is similar to the function of the same
  !              name from the tensbs library without the work array
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, parametric position 1
  !     v       - Real, parametric position 2
  !     w       - Real, parametric position 3
  !     idu     - Integer, derivative in u to evaluate
  !     idv     - Integer, derivative in v to evaluate
  !     idw     - Integer, derivative in w to evaluate
  !     tu      - Real, U Knot vector. Length nctlu+ku
  !     tv      - Real, V Knot vector. Length nctlv+kv
  !     tv      - Real, W Knot vector. Length nctlw+kw
  !     Nctlu   - Integer, Number of u control points
  !     Nctlv   - Integer, Number of v control points
  !     Nctlw   - Integer, Number of w control points
  !     ku      - Integer, order of spline in u direction
  !     kv      - Integer, order of spline in v direction
  !     kw      - Integer, order of spline in w direction
  !     coef    - Complex, Array size Nctlu x Nctlv x Nctlw, B-spline coefficients
  !
  !     Ouput 
  !     b3val   - Complex, Evaluated point or derivative

  use precision
  implicit none

  ! Input
  integer         , intent(in)        :: ku, kv, kw, nctlu, nctlv, nctlw
  integer         , intent(in)        :: idu, idv, idw
  real(kind=realType), intent(in)     :: u, v, w
  real(kind=realType), intent(in)     :: tu(nctlu+ku), tv(nctlv+kv), tw(nctlw+kw)
  complex(kind=realType), intent(in)  :: coef(nctlu,nctlv,nctlw)

  ! Output
  complex(kind=realType)           :: cb3val

  ! Working
  integer                          :: i, j, k
  integer                          :: ileftu, ileftv, ileftw
  integer                          :: istartu, istartv, istartw
  real(kind=realType)              :: Bu(ku), Bv(kv), Bw(kw)
  real(kind=realType)              :: Bud(ku, ku), Bvd(kv, kv), Bwd(kw, kw)

  ! Find knot spans
  call findSpan(u, ku, tu, nctlu, ileftu)
  istartu = ileftu-ku

  call findSpan(v, kv, tv, nctlv, ileftv)
  istartv = ileftv-kv

  call findSpan(w, kw, tw, nctlw, ileftw)
  istartw = ileftw-kw

  cb3val = cmplx(0.0, 0.0)

  ! Check if we can just use regular basis functions:
  if (idu + idv + idw == 0) then
     call basis(tu, nctlu, ku, u, ileftu, Bu)
     call basis(tv, nctlv, kv, v, ileftv, Bv)
     call basis(tw, nctlw, kw, w, ileftw, Bw)

     do i=1, ku
        do j=1, kv
           do k=1, kw
              cb3val = cb3val + Bu(i)*Bv(j)*Bw(k)*&
                   coef(istartu+i, istartv+j, istartw+k)
           end do
        end do
     end do
  else
     ! We need derivative basis functions
     call derivBasis(tu, nctlu, ku, u, ileftu, idu, Bu)
     call derivBasis(tv, nctlv, kv, v, ileftv, idv, Bv)
     call derivBasis(tw, nctlw, kw, w, ileftw, idw, Bw)

     do i=1, ku
        do j=1, kv
           do k=1, kw
              cb3val = cb3val + Bud(idu+1, i)*Bvd(idv+1, j)*Bwd(idw+1, k) * &
                   coef(istartu+i, istartv+j, istartw+k)
           end do
        end do
     end do
  end if

end function cb3val
