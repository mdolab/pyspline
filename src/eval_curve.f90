subroutine eval_curve(s, t, k, coef, nctl, ndim, n, val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surf is a vector version of the B-spline evaluation function
  !
  !     Description of Arguments
  !     Input
  !     s       - Real, Vector of s coordinates, length n
  !     t       - Real, Knot vector. Length nctl+k
  !     k       - Integer, order of B-spline 
  !     coef    - Real, Array of b-sline coefficients  Size (ndim, nctl)
  !     nctl    - Integer, Number of control points
  !
  !     Ouput 
  !     val     - Real, Evaluated points, size ndim by n

  use precision
  implicit none
  ! Input
  integer         , intent(in)     :: k, nctl, ndim, n
  real(kind=realType), intent(in)  :: s(n)
  real(kind=realType), intent(in)  :: t(nctl+k)
  real(kind=realType), intent(in)  :: coef(ndim, nctl)

  ! Output
  real(kind=realType), intent(out) :: val(ndim, n)

  ! Working
  integer                          :: i, l, idim, istart, ileft
  real(kind=realType)              :: basisu(k)

  val(:, :) = 0.0
  do i=1, n
     call findSpan(s(i), k, t, nctl, ileft)
     call basis(t, nctl, k, s(i), ileft, basisu)
     istart = ileft-k
     do l=1, k
        do idim=1, ndim
           val(idim, i) = val(idim, i) + basisu(l)*coef(idim, istart+l)
        end do
     end do
  end do
end subroutine eval_curve

subroutine eval_curve_deriv(s, t, k, coef, nctl, ndim, val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surf_deriv is a scalar version of the 
  !              B-spline derivative evaluation function
  !
  !     Description of Arguments
  !     Input
  !     s       - Real, s coordinate
  !     t       - Real, Knot vector. Length nctl+k
  !     k       - Integer, order of B-spline
  !     coef    - Real, Array of B-spline coefficients. Size (ndim, nctl)
  !     nctl   - Integer, Number of control points
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim

  use precision
  implicit none
 
  ! Input
  integer         , intent(in)          :: k, nctl, ndim
  real(kind=realType), intent(in)       :: s
  real(kind=realType), intent(in)       :: t(nctl+k)
  real(kind=realType), intent(in)       :: coef(ndim, nctl)

  ! Output
  real(kind=realType), intent(out)      :: val(ndim)

  ! Working
  integer                               :: idim

  ! Functions
  real(kind=realType)                   :: bvalu

  do idim=1, ndim
     val(idim) = bvalu(t, coef(idim, :), nctl, k, 1, s)
  end do

end subroutine eval_curve_deriv

subroutine eval_curve_deriv2(s, t, k, coef, nctl, ndim, val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_curve_deriv2 evals the 2nd derivative
  !              B-spline derivative evaluation function
  !
  !     Description of Arguments
  !     Input
  !     s       - Real, s coordinate
  !     t       - Real, Knot vector. Length nctl+k
  !     k       - Integer, order of B-spline
  !     coef    - Real, Array of B-spline coefficients Size (ndim, nctl)
  !     nctl    - Integer, Number of control points
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim

  use precision
  implicit none
  ! Input
  integer         , intent(in)             :: k, nctl, ndim
  real(kind=realType), intent(in)          :: s
  real(kind=realType), intent(in)          :: t(nctl+k)
  real(kind=realType), intent(in)          :: coef(ndim, nctl)

  ! Output
  real(kind=realType), intent(out)         :: val(ndim)

  ! Working
  integer                                  :: idim

  ! Functions
  real(kind=realType)                      :: bvalu

  if (k == 2) then
     val(:) = 0.0
  else
     do idim=1, ndim
        val(idim) = bvalu(t, coef(idim, :), nctl, k, 2, s)
     end do
  end if

end subroutine eval_curve_deriv2


subroutine eval_curve_c(s, t, k, coef, nctl, ndim, n, val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_curve_complex is a special version for evaluation
  !     of data when control points are used in a complex-step analysis
  !     The actual derivative is computed analytically, but the output looks
  !     like the CS went all the way through
  !
  !     Description of Arguments
  !     Input
  !     s       - Real, s coordinate
  !     t       - Real, Knot vector. Length nctl+k
  !     k       - Integer, order of B-spline
  !     coef    - Real, Array of B-spline coefficients  Size (ndim, nctl)
  !     nctl   - Integer, Number of control points
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim

  use precision
  implicit none
  ! Input
  integer            , intent(in)       :: k, nctl, ndim, n
  real(kind=realType), intent(in)       :: s(n)
  real(kind=realType), intent(in)       :: t(nctl+k)
  complex*16         , intent(in)       :: coef(ndim, nctl)

  ! Output
  complex*16      , intent(out)         :: val(ndim, n)

  ! Working
  integer                               :: i, l, istart, ileft, idim
  real(kind=realType)                   :: basisu(k)

  val(:, :) = 0.0
  do i=1, n
     call findSpan(s(i), k, t, nctl, ileft)
     call basis(t, nctl, k, s(i), ileft, basisu)

     istart = ileft - k
     do l=1, k
        do idim=1, ndim
           val(idim, i) = val(idim, i) + &
                cmplx(basisu(l)*real(coef(idim, istart+l)), &
                basisu(l)*aimag(coef(idim, istart+l)))

        end do
     end do
  end do
end subroutine eval_curve_c

subroutine eval_curve_deriv_c(s, t, k, coef, nctl, ndim, val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surf_deriv is a scalar version of the 
  !              B-spline derivative evaluation function
  !
  !     Description of Arguments
  !     Input
  !     s       - Real, s coordinate
  !     t       - Real, Knot vector. Length nctl+k
  !     k       - Integer, order of B-spline
  !     coef    - Real, Array of B-spline coefficients. Size (ndim, nctl)
  !     nctl   - Integer, Number of control points
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim

  use precision
  implicit none
  ! Input
  integer   , intent(in)          :: k, nctl, ndim
  complex*16, intent(in)          :: s
  complex*16, intent(in)          :: t(nctl+k)
  complex*16, intent(in)          :: coef(ndim, nctl)

  ! Output
  complex*16, intent(out)         :: val(ndim)

  ! Working
  integer                         :: idim, inbv
  complex*16                      :: work(3*k)

  ! Functions
  complex*16                      :: bvalu_c

  inbv = 1

!   do idim=1, ndim
!      val(idim) = bvalu_c(t, coef(idim, :), nctl, k, 1, s, inbv, work)
!   end do

end subroutine eval_curve_deriv_c

subroutine eval_curve_deriv2_c(s, t, k, coef, nctl, ndim, val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_curve_deriv2 evals the 2nd derivative
  !              B-spline derivative evaluation function
  !
  !     Description of Arguments
  !     Input
  !     s       - Real, s coordinate
  !     t       - Real, Knot vector. Length nctl+k
  !     k       - Integer, order of B-spline
  !     coef    - Real, Array of B-spline coefficients Size (ndim, nctl)
  !     nctl    - Integer, Number of control points
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim

  use precision
  implicit none
  ! Input
  integer   , intent(in)          :: k, nctl, ndim
  complex*16, intent(in)          :: s
  complex*16, intent(in)          :: t(nctl+k)
  complex*16, intent(in)          :: coef(ndim, nctl)

  ! Output
  complex*16, intent(out)         :: val(ndim)

  ! Working
  integer                         :: idim, inbv
  complex*16                      :: work(3*k)

  ! Functions
  complex*16                      :: bvalu

  inbv = 1

  if (k == 2) then
     val(:) = 0.0
  else
     do idim=1, ndim
        val(idim) = bvalu(t, coef(idim, :), nctl, k, 2, s)
     end do
  end if

end subroutine eval_curve_deriv2_c
