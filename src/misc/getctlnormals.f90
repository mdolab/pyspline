subroutine getctlnormals(coef,tx,ty,kx,ky,nx,ny,normals)
!***DESCRIPTION
!
!     Written by Gaetan Kenway
!
!     Abstract: getctlnormals produces the surface normals for each
!         control point in coef. The idea is to use these normals to
!         move control points
!
!     Description of Arguments
!     Input
!     coef    - Real, Matrix of coefficients: size: nx by ny by 3
!     tx      - Real, Vector: Knots for x, size nx + k
!     ty      - Real, Vector: Knots for y, size ny + k
!     kx      - Integer: Order of spline in x
!     ky      - Integer: Order of spline in y
!     nx      - Integer: Number of coefficients in x
!     ny      - Integer: Number of coefficients in y
!
!     Ouput
!     normals - Real, Matrix of normals: size: nx by ny by 3
!
  implicit none

  integer nx,ny,i,j,idim
  integer,          intent(in)     :: kx,ky
  double precision, intent(in)     :: tx(kx+nx)
  double precision, intent(in)     :: ty(kx+ny)
  double precision, intent(inout)  :: coef(nx,ny,3)
  double precision, intent(out)    :: normals(nx,ny,3)

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
           dx(idim) = b2val(tx(kx-1 + i),ty(ky-1 + j),1,0,tx,ty,nx,ny,kx,ky,coef(:,:,idim),work)
           dy(idim) = b2val(tx(kx-1 + i),ty(ky-1 + j),0,1,tx,ty,nx,ny,kx,ky,coef(:,:,idim),work)
        end do

        ! Now get the cross product 

        n(1) = dx(2)*dy(3)-dx(3)*dy(2)
        n(2) = dx(3)*dy(1)-dx(1)*dy(3)
        n(3) = dx(1)*dy(2)-dx(2)*dy(1)

        ! Normalize  the cross product 
        n= n /sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))

        normals(i,j,:) = n
     end do
  end do
end subroutine getctlnormals



     
