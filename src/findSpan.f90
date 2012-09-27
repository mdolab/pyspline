subroutine findSpan(u, k, t, nctl, ind)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway. Adapted from "The NURBS BooK" Algorithm A2.1
  !
  !     Abstract: Determine the knot span index

  !     Description of Arguments
  !     Input:
  !     u       - Real, parametric location we are looking for
  !     k       - Integer, order of B-spline 
  !     t       - Real, size(nctl + k), knot vector
  !
  !     Ouput: 
  !     ind     - Integer, knot span index

  use precision
  implicit none

  ! Input 
  integer, intent(in)             :: k, nctl
  real(kind=realType), intent(in) :: u, t(nctl+k)
  ! Output
  integer, intent(out)            :: ind

  ! Working
  integer :: low, mid, high

  if (u == t(nctl+1)) then
     ind = nctl
  else
     low = k
     high = nctl+1

     ! Do a binary search
     mid = (low+high)/2
     do while ( (u < t(mid) .or. u >= t(mid+1)))
        if (u < t(mid)) then
           high = mid
        else
           low = mid
        end if

        mid = (low+high)/2
     end do

     ind = mid
  end if

end subroutine findSpan
