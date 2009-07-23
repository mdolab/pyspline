subroutine updateSurfacePoints(coef,update,tx,ty,kx,ky,nx,ny)


!*** DESCRIPTION
!
!     Written by Gaetan Kenway
! 
!     The function takes the defination of a bivariate spline surface
!     as defined by coef,tx,ty,kx and ky and displaces the control
!     points by the values contained in update along the vector normal
!     to the surface. 
!
!     Description of Arguments:
!     Input:
!     coef    - Real, Matrix of coefficients, size nx by ny by 3 
!               (This only works with 3 spatial degrees of freedom)
!     update  - Real, Matrix of displacements, size, nx by ny
!     tx      - Real, vector of knots in x, size nx + kx
!     ty      - Real, vector of knots in y, size ny + ky
!     nx      - Integer, The number of control points in x
!     ny      - Integer, The number of control points in y
!
!     Output:
!     coef    - Real, Matrix of displaced coefficients, size nx by ny
!
  implicit none

  integer nx,ny,i,j,idim
  integer,          intent(in)     :: kx,ky
  double precision, intent(in)     :: tx(kx+nx)
  double precision, intent(in)     :: ty(kx+ny)
  double precision, intent(inout)  :: coef(nx,ny,3)
  double precision, intent(in)     :: update(nx,ny)

  double precision                 :: dx(3)
  double precision                 :: dy(3)
  double precision                 :: n(3)
  double precision                 :: WORK(3*MAX(KX,KY)+KY)

  double precision b2val

  do i =1,Nx
     do j = 1,Ny

        ! Get the normal for the current point:
        
        ! First get the directional derivatives 
        do idim = 1,3
           dx(idim) = b2val(tx(kx-2 + i),ty(ky-2 + j),1,0,tx,ty,nx,ny,kx,ky,coef(:,:,idim),work)
           dy(idim) = b2val(tx(kx-2 + i),ty(ky-2 + j),0,1,tx,ty,nx,ny,kx,ky,coef(:,:,idim),work)
        end do

        ! Now get the cross product 

        n(1) = dx(2)*dy(3)-dx(3)*dy(2)
        n(2) = dx(3)*dy(1)-dx(1)*dy(3)
        n(3) = dx(1)*dy(2)-dx(2)*dy(1)

        ! Normalizet the cross product
        n= n /sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))

        ! Move control point along the normal vector
        do idim = 1,3
           coef(i,j,idim) = coef(i,j,idim) - update(i,j)*n(idim)
        end do
     end do
  end do
end subroutine updateSurfacePoints
