subroutine tfi2d(e0, e1, e2, e3, Nu, Nv, X)
 
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: Perform a simple 2D transfinite interpolation in 3
  !               spatial dimensions. This function is not used
  !               directly in pySpline, but from pyGeo. 
  !
  !     Description of Arguments
  !     Input
  !     e0      - Real, Vector or coordinates along 0th edge. Size (3, Nu)
  !     e1      - Real, Vector or coordinates along 1st edge. Size (3, Nu)
  !     e2      - Real, Vector or coordinates along 2nd edge. Size (3, Nv)
  !     e3      - Real, Vector or coordinates along 3rd edge. Size (3, Nv)
  !
  !     Ouput 
  !     X       - Real, Evaluated points, size 3 x Nu x Nv
  ! 
  
  implicit none

  integer, intent(in) :: Nu,Nv
  double precision, intent(in) :: e0(3,Nu)
  double precision, intent(in) :: e1(3,Nu)

  double precision, intent(in) :: e2(3,Nv)
  double precision, intent(in) :: e3(3,Nv)

  double precision, intent(out) :: X(3,Nv,Nu)

  double precision :: U(Nu),V(Nv)
  integer :: i,j

  do i=1,Nu
     U(i) = dble(i-1)/(Nu-1)
  end do

  do j=1,Nv
     V(j) = dble(j-1)/(Nv-1)
  end do

  do i=1,Nu
     do j=1,Nv
        X(:,j,i) = (1-V(j))*e0(:,i) + V(j)*e1(:,i) + &
                   (1-U(i))*e2(:,j) + U(i)*e3(:,j) - ( &
                   (  U(i))*(  V(j))*e1(:,Nu) + &
                   (  U(i))*(1-V(j))*e0(:,Nu) + &
                   (1-U(i))*(  V(j))*e1(:,1)  + &
                   (1-U(i))*(1-V(j))*e0(:,1) &
                   )
     end do
  end do

end subroutine tfi2d
