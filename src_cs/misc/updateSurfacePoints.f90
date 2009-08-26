subroutine updateSurfacePoints(coef,l_index,g_index,update,s_coef,tx,ty,kx,ky,nx,ny,ncoef,n)


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
!     coef    - Real, List of global coefficients, size n by 3
!     l_index - Integer, List of x-y positions on surface each control point 
!               refers to. size, n by 2
!     g_index - Integer, List of positions of 'coef' in the global coef array
!     s_coef  - Real, matrix size nx by ny: Coefficients for the spline surface
!     tx      - Real, vector of knots in x, size nx + kx
!     ty      - Real, vector of knots in y, size ny + ky
!     nx      - Integer, The number of control points in x
!     ny      - Integer, The number of control points in y
!
!     Output:
!     coef    - Updated global coefficient list  
!
  implicit none

  integer nx,ny,i,j,idim,icoef,n,ncoef
  ! Input
  complex*16, intent(inout)  :: coef(ncoef,3)
  integer   , intent(in)     :: l_index(n,2)
  integer   , intent(in)     :: g_index(n)
  complex*16, intent(in)     :: update(n)
  complex*16, intent(in)     :: s_coef(nx,ny,3)
  integer,    intent(in)     :: kx,ky
  complex*16, intent(in)     :: tx(kx+nx)
  complex*16, intent(in)     :: ty(ky+ny)

  ! Local
  complex*16                 :: dx(3)
  complex*16                 :: dy(3)
  complex*16                 :: normal(3)
  complex*16                 :: WORK(3*MAX(KX,KY)+KY)

  complex*16 b2val

  do icoef = 1,N
     
     ! Get the normal for the current point:
     ! First get the directional derivatives 
     i = l_index(icoef,1) + 1 !plus 1 for python->fortran conversion
     j = l_index(icoef,2) + 1 !plus 1 for python->fortran conversion

     do idim = 1,3
        dx(idim) = b2val(tx(kx-1 + i),ty(ky-1 + j),1,0,tx,ty,nx,ny,kx,ky,s_coef(:,:,idim),work)
        dy(idim) = b2val(tx(kx-1 + i),ty(ky-1 + j),0,1,tx,ty,nx,ny,kx,ky,s_coef(:,:,idim),work)
     end do

     ! Now get the cross product 

     normal(1) = dx(2)*dy(3)-dx(3)*dy(2)
     normal(2) = dx(3)*dy(1)-dx(1)*dy(3)
     normal(3) = dx(1)*dy(2)-dx(2)*dy(1)

     ! Normalizet the cross product
     normal=normal/sqrt(normal(1)*normal(1)+normal(2)*normal(2)+normal(3)*normal(3))

     ! Move control point along the normal vector
     do idim = 1,3
        coef(g_index(icoef+1),idim) = coef(g_index(icoef+1),idim) - update(icoef)*normal(idim)
     end do
  end do
end subroutine updateSurfacePoints
