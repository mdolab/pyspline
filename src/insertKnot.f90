subroutine insertKnot(u, r, t, k, coef, nctl, ndim, t_new, coef_new, ileft)

    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Abstract insertKnot inserts a knot u into the curve, r times
    !     Adapted from "The NURBS Book" Algorithm 5.1
    !     Description of Arguments
    !     Input
    !     u       - Real, location of knot to insert
    !     r       - Integer, Insert r times
    !     t       - Real,Knot vector. Length nctl+k
    !     k       - Integer,order of B-spline
    !     coef    - Real,Array of B-spline coefficients  Size (ndim,nctl)
    !     nctl    - Integer,Number of control points
    !     ndim    - Integer, dimension of curve

    !     Ouput
    !     t_new    - Real, vector of lenght(nctl+k+r)
    !     coef_new - Real, Array of new cofficients size(ndim,nctl+r)
    !     ileft    - Integer of position of knot insertion

    use precision
    implicit none

    ! Input
    integer, intent(inout) :: r
    integer, intent(in) :: k, nctl, ndim
    real(kind=realType), intent(in) :: u
    real(kind=realType), intent(in) :: t(nctl + k)
    real(kind=realType), intent(in) :: coef(ndim, nctl)

    ! Output
    real(kind=realType), intent(out) :: t_new(nctl + k + r)
    real(kind=realType), intent(out) :: coef_new(ndim, nctl + r)
    integer, intent(out) :: ileft

    ! Working
    integer :: s, i, j, L
    real(kind=realType) :: alpha, temp(ndim, k)

    call findSpan(u, k, t, nctl, ileft)

    ! Compute its multiplicity
    s = 0 ! Knot multiplicity
    mult: do i = 0, k - 1
        if (abs(t(ileft - i) - u) < 1e-12) then
            s = s + 1
        else
            exit mult
        end if
    end do mult

    ! We need to make sure that the requested multipliity r, plus
    ! the actual multiplicity of this know is less than (k-1)

    if (s + r + 1 > k) then
        r = k - 1 - s
    end if

    ! If we *can't* insert the knot, we MUST copy t and coef to t_new
    ! and coef_new and return
    if (r == 0) then
        coef_new(:, 1:nctl) = coef(:, 1:nctl)
        t_new(1:nctl + k) = t(1:nctl + k)
        return
    end if

    ! --------- Load New Knot Vector -------
    do i = 1, ileft
        t_new(i) = t(i)
    end do

    do i = 1, r
        t_new(ileft + i) = u
    end do

    do i = ileft + 1, nctl + k
        t_new(i + r) = t(i)
    end do

    ! -------- Save unaltered Control Points
    do i = 1, ileft - (k - 1)
        coef_new(:, i) = coef(:, i)
    end do

    do i = ileft - s, nctl
        coef_new(:, i + r) = coef(:, i)
    end do

    do i = 0, k - 1 - s
        temp(:, i + 1) = coef(:, ileft - (k - 1) + i)
    end do

    do j = 1, r
        L = ileft - (k - 1) + j
        do i = 0, k - 1 - j - s
            alpha = (u - t(L + i)) / (t(i + ileft + 1) - t(L + i))
            temp(:, i + 1) = alpha * temp(:, i + 2) + (1.0 - alpha) * temp(:, i + 1)
        end do
        coef_new(:, L) = temp(:, 1)
        coef_new(:, ileft + r - j - s) = temp(:, k - j - s)
    end do

    do i = L + 1, ileft - s - 1
        coef_new(:, i) = temp(:, i - L + 1)
    end do

    call findSpan(u, k, t_new, nctl + r, ileft)
end subroutine insertKnot

