subroutine eval_surface(u, v, tu, tv, ku, kv, coef, nctlu, nctlv, ndim, &
                        n, m, val)

    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Abstract eval_surface evaluates (possibly) many points on the surface
    !
    !     Description of Arguments
    !     Input
    !     u       - Real, u coordinate, size(m, n)
    !     v       - Real, v coordinate, size(m, n)
    !     tu      - Real, Knot vector in u. size(nctlu+ku)
    !     tv      - Real, Knot vector in v. size(nctlv+kv)
    !     ku      - Integer, order of B-spline in u
    !     kv      - Integer, order of B-spline in v
    !     coef    - Real, Array of B-spline coefficients  Size (ndim, nctlv, nctlu)
    !     nctlu   - Integer, Number of control points in u
    !     nctlv   - Integer, Number of control points in v
    !     ndim    - Integer, Spatial Dimension
    !
    !     Ouput
    !     val     - Real, Evaluated point(s), size (ndim, m, n)

    use precision
    implicit none

    ! Input
    integer, intent(in) :: ku, kv, nctlu, nctlv, ndim, n, m
    real(kind=realType), intent(in) :: u(m, n), v(m, n)
    real(kind=realType), intent(in) :: tu(nctlu + ku), tv(nctlv + kv)
    real(kind=realType), intent(in) :: coef(ndim, nctlv, nctlu)

    ! Output
    real(kind=realType), intent(out) :: val(ndim, m, n)

    ! Working
    integer :: idim, istartu, istartv, i, j, ii, jj
    integer :: ileftu, ileftv
    real(kind=realType) :: basisu(ku), basisv(kv)

    val(:, :, :) = 0.0
    do ii = 1, n
        do jj = 1, m
            ! U
            call findSpan(u(jj, ii), ku, tu, nctlu, ileftu)
            call basis(tu, nctlu, ku, u(jj, ii), ileftu, basisu)
            istartu = ileftu - ku

            ! V
            call findSpan(v(jj, ii), kv, tv, nctlv, ileftv)
            call basis(tv, nctlv, kv, v(jj, ii), ileftv, basisv)
            istartv = ileftv - kv

            do i = 1, ku
                do j = 1, kv
                    do idim = 1, ndim
                    val(idim, jj, ii) = val(idim, jj, ii) + basisu(i) * basisv(j) * coef(idim, istartv + j, istartu + i)
                    end do
                end do
            end do
        end do
    end do

end subroutine eval_surface

subroutine eval_surface_deriv(u, v, tu, tv, ku, kv, coef, nctlu, nctlv, &
                              ndim, val)

    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Abstract eval_surface_deriv evaluates the derivative of a
    !     point on the surface
    !
    !     Description of Arguments
    !     Input
    !     u       - Real, u coordinate
    !     v       - Real, v coordinate
    !     tu      - Real, Knot vector in u. size(nctlu+ku)
    !     tv      - Real, Knot vector in v. size(nctlv+kv)
    !     ku      - Integer, order of B-spline in u
    !     kv      - Integer, order of B-spline in v
    !     coef    - Real, Array of B-spline coefficients  Size (ndim, nctlv, nctlu)
    !     nctlu   - Integer, Number of control points in u
    !     nctlv   - Integer, Number of control points in v
    !     ndim    - Integer, Spatial Dimension
    !
    !     Ouput
    !     val     - Real, Evaluated derivatives, size (ndim, 2)

    use precision
    implicit none

    ! Input
    integer, intent(in) :: ku, kv, nctlu, nctlv, ndim
    real(kind=realType), intent(in) :: u, v
    real(kind=realType), intent(in) :: tu(nctlu + ku), tv(nctlv + kv)
    real(kind=realType), intent(in) :: coef(ndim, nctlv, nctlu)

    ! Output
    real(kind=realType), intent(out) :: val(ndim, 2)

    ! Working
    integer :: idim
    real(kind=realType) :: b2val

    do idim = 1, ndim
        val(idim, 1) = b2val(v, u, 0, 1, tv, tu, nctlv, nctlu, kv, ku, coef(idim, :, :))
        val(idim, 2) = b2val(v, u, 1, 0, tv, tu, nctlv, nctlu, kv, ku, coef(idim, :, :))
    end do

end subroutine eval_surface_deriv

subroutine eval_surface_deriv2(u, v, tu, tv, ku, kv, coef, nctlu, nctlv, ndim, val)

    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Abstract eval_surface_deriv2 evaluates the second derivative
    !     of a point on the surface
    !
    !     Description of Arguments
    !     Input
    !     u       - Real, u coordinate
    !     v       - Real, v coordinate
    !     tu      - Real, Knot vector in u. size(nctlu+ku)
    !     tv      - Real, Knot vector in v. size(nctlv+kv)
    !     ku      - Integer, order of B-spline in u
    !     kv      - Integer, order of B-spline in v
    !     coef    - Real, Array of B-spline coefficients  Size (ndim, nctlv, nctlu)
    !     nctlu   - Integer, Number of control points in u
    !     nctlv   - Integer, Number of control points in v
    !     ndim    - Integer, Spatial Dimension
    !
    !     Ouput
    !     val     - Real, Evaluated second derivatives, size (ndim, 2, 2)

    use precision
    implicit none

    ! Input
    integer, intent(in) :: ku, kv, nctlu, nctlv, ndim
    real(kind=realType), intent(in) :: u, v
    real(kind=realType), intent(in) :: tu(nctlu + ku), tv(nctlv + kv)
    real(kind=realType), intent(in) :: coef(ndim, nctlv, nctlu)

    ! Output
    real(kind=realType), intent(out) :: val(ndim, 2, 2)

    ! Working
    integer :: idim
    real(kind=realType) :: b2val

    do idim = 1, ndim
        val(idim, 1, 1) = b2val(v, u, 0, 2, tv, tu, nctlv, nctlu, kv, ku, coef(idim, :, :))
        val(idim, 1, 2) = b2val(v, u, 1, 1, tv, tu, nctlv, nctlu, kv, ku, coef(idim, :, :))
        val(idim, 2, 1) = val(idim, 1, 2)
        val(idim, 2, 2) = b2val(v, u, 2, 0, tv, tu, nctlv, nctlu, kv, ku, coef(idim, :, :))
    end do

end subroutine eval_surface_deriv2
