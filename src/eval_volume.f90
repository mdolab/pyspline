subroutine eval_volume(u, v, w, tu, tv, tw, ku, kv, kw, coef, &
                       nctlu, nctlv, nctlw, ndim, n, m, l, val)

    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Abstract eval_volume evaluates (possibly) many points on the
    !              b-spline volume
    !
    !     Description of Arguments
    !     Input
    !     u       - Real, u coordinate, size(l, m, n)
    !     v       - Real, v coordinate, size(l, m, n)
    !     w       - Real, w coordinate, size(l, m, n)
    !     tu      - Real, Knot vector in u. Length nctlu+ku
    !     tv      - Real, Knot vector in v. Length nctlv+kv
    !     tw      - Real, Knot vector in w. Length nctlv+kw
    !     ku      - Integer, order of B-spline in u
    !     kv      - Integer, order of B-spline in v
    !     kw      - Integer, order of B-spline in w
    !     coef    - Real, Array of B-spline coefficients
    !                 Size (ndim, nctlw, nctlv, nctlu)
    !     nctlu   - Integer, Number of control points in u
    !     nctlv   - Integer, Number of control points in v
    !     nctlw   - Integer, Number of control points in w
    !     ndim    - Integer, Spatial Dimension
    !
    !     Ouput
    !     val     - Real, Evaluated points, size (ndim, l, m, n)

    use precision
    implicit none

    ! Input
    integer, intent(in) :: ku, kv, kw, nctlu, nctlv, nctlw
    integer, intent(in) :: ndim, n, m, l
    real(kind=realType), intent(in) :: u(l, m, n), v(l, m, n), w(l, m, n)
    real(kind=realType), intent(in) :: tu(nctlu + ku), tv(nctlv + kv), tw(nctlw + kw)
    real(kind=realType), intent(in) :: coef(ndim, nctlw, nctlv, nctlu)

    ! Output
    real(kind=realType), intent(out) :: val(ndim, l, m, n)

    ! Working
    integer :: idim, istartu, istartv, istartw
    integer :: i, j, k, ii, jj, kk
    integer :: ileftu, ileftv, ileftw
    real(kind=realType) :: basisu(ku), basisv(kv), basisw(kw)

    val(:, :, :, :) = 0.0
    do ii = 1, n
        do jj = 1, m
            do kk = 1, l
                ! U
                call findSpan(u(kk, jj, ii), ku, tu, nctlu, ileftu)
                call basis(tu, nctlu, ku, u(kk, jj, ii), ileftu, basisu)
                istartu = ileftu - ku

                ! V
                call findSpan(v(kk, jj, ii), kv, tv, nctlv, ileftv)
                call basis(tv, nctlv, kv, v(kk, jj, ii), ileftv, basisv)
                istartv = ileftv - kv

                ! W
                call findSpan(w(kk, jj, ii), kw, tw, nctlw, ileftw)
                call basis(tw, nctlw, kw, w(kk, jj, ii), ileftw, basisw)
                istartw = ileftw - kw

                do i = 1, ku
                    do j = 1, kv
                        do k = 1, kw
                            do idim = 1, ndim
                                val(idim, kk, jj, ii) = val(idim, kk, jj, ii) + &
                                                        basisu(i) * basisv(j) * basisw(k) * &
                                                        coef(idim, istartw + k, istartv + j, istartu + i)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end subroutine eval_volume

subroutine eval_volume_deriv(u, v, w, tu, tv, tw, ku, kv, kw, coef, &
                             nctlu, nctlv, nctlw, ndim, val)

    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Abstract eval_volume_deriv evaluates the first derivative on a
    !              volume
    !
    !     Description of Arguments
    !     Input
    !     u       - Real, u coordinate
    !     v       - Real, v coordinate
    !     w       - Real, w coordinate
    !     tu      - Real, Knot vector in u. Length nctlu+ku
    !     tv      - Real, Knot vector in v. Length nctlv+kv
    !     tw      - Real, Knot vector in w. Length nctlv+kw
    !     ku      - Integer, order of B-spline in u
    !     kv      - Integer, order of B-spline in v
    !     kw      - Integer, order of B-spline in w
    !     coef    - Real, Array of B-spline coefficients
    !                 Size (ndim, nctlw, nctlv, nctlu)
    !     nctlu   - Integer, Number of control points in u
    !     nctlv   - Integer, Number of control points in v
    !     nctlw   - Integer, Number of control points in w
    !     ndim    - Integer, Spatial Dimension
    !
    !     Ouput
    !     val     - Real, Evaluated derivatvie, size(ndim, 3)

    use precision
    implicit none

    ! Input
    integer, intent(in) :: ku, kv, kw, nctlu, nctlv, nctlw, ndim
    real(kind=realType), intent(in) :: u, v, w
    real(kind=realType), intent(in) :: tu(nctlu + ku), tv(nctlv + kv), tw(nctlw + kw)
    real(kind=realType), intent(in) :: coef(ndim, nctlw, nctlv, nctlu)

    ! Output
    real(kind=realType), intent(out) :: val(ndim, 3)

    ! Working
    integer :: idim
    real(kind=realType) :: b3val

    do idim = 1, ndim
        val(idim, 1) = b3val(w, v, u, 0, 0, 1, tw, tv, tu, nctlw, nctlv, nctlu, &
                             kw, kv, ku, coef(idim, :, :, :))

        val(idim, 2) = b3val(w, v, u, 0, 1, 0, tw, tv, tu, nctlw, nctlv, nctlu, &
                             kw, kv, ku, coef(idim, :, :, :))

        val(idim, 3) = b3val(w, v, u, 1, 0, 0, tw, tv, tu, nctlw, nctlv, nctlu, &
                             kw, kv, ku, coef(idim, :, :, :))
    end do
end subroutine eval_volume_deriv

subroutine eval_volume_deriv2(u, v, w, tu, tv, tw, ku, kv, kw, coef, &
                              nctlu, nctlv, nctlw, ndim, val)

    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Abstract eval_volume_deriv2 evaluates the second derivative on a
    !              volume
    !
    !     Description of Arguments
    !     Input
    !     u       - Real, u coordinate
    !     v       - Real, v coordinate
    !     w       - Real, w coordinate
    !     tu      - Real, Knot vector in u. Length nctlu+ku
    !     tv      - Real, Knot vector in v. Length nctlv+kv
    !     tw      - Real, Knot vector in w. Length nctlv+kw
    !     ku      - Integer, order of B-spline in u
    !     kv      - Integer, order of B-spline in v
    !     kw      - Integer, order of B-spline in w
    !     coef    - Real, Array of B-spline coefficients
    !                 Size (ndim, nctlw, nctlv, nctlu)
    !     nctlu   - Integer, Number of control points in u
    !     nctlv   - Integer, Number of control points in v
    !     nctlw   - Integer, Number of control points in w
    !     ndim    - Integer, Spatial Dimension
    !
    !     Ouput
    !     val     - Real, Evaluated derivatvie, size(ndim, 3, 3)

    use precision
    implicit none
    ! Input
    integer, intent(in) :: ku, kv, kw, nctlu, nctlv, nctlw, ndim
    real(kind=realType), intent(in) :: u, v, w
    real(kind=realType), intent(in) :: tu(nctlu + ku), tv(nctlv + kv), tw(nctlw + kw)
    real(kind=realType), intent(in) :: coef(ndim, nctlw, nctlv, nctlu)

    ! Output
    real(kind=realType), intent(out) :: val(ndim, 3, 3)

    ! Working
    integer :: idim

    real(kind=realType) :: b3val

    do idim = 1, ndim

        ! Row 1
        if (ku >= 3) then
            val(idim, 1, 1) = b3val(w, v, u, 0, 0, 2, tw, tv, tu, nctlw, nctlv, &
                                    nctlu, kw, kv, ku, coef(idim, :, :, :))
        else
            val(idim, 1, 1) = 0.0
        end if
        val(idim, 1, 2) = b3val(w, v, u, 0, 1, 1, tw, tv, tu, nctlw, nctlv, nctlu, &
                                kw, kv, ku, coef(idim, :, :, :))
        val(idim, 1, 3) = b3val(w, v, u, 1, 0, 1, tw, tv, tu, nctlw, nctlv, nctlu, &
                                kw, kv, ku, coef(idim, :, :, :))

        ! Row 2
        val(idim, 2, 1) = val(idim, 1, 2)
        if (kv >= 3) then
            val(idim, 2, 2) = b3val(w, v, u, 0, 2, 0, tw, tv, tu, nctlw, nctlv, &
                                    nctlu, kw, kv, ku, coef(idim, :, :, :))
        else
            val(idim, 2, 2) = 0.0
        end if
        val(idim, 2, 3) = b3val(w, v, u, 1, 1, 0, tw, tv, tu, nctlw, nctlv, nctlu, &
                                kw, kv, ku, coef(idim, :, :, :))

        ! Row 3
        val(idim, 3, 1) = val(idim, 1, 3)
        val(idim, 3, 2) = val(idim, 2, 3)
        if (kw >= 3) then
            val(idim, 3, 3) = b3val(w, v, u, 2, 0, 0, tw, tv, tu, nctlw, nctlv, &
                                    nctlu, kw, kv, ku, coef(idim, :, :, :))
        else
            val(idim, 3, 3) = 0.0
        end if
    end do

end subroutine eval_volume_deriv2
