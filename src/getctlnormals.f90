subroutine getctlnormals(coef,l_index,g_index,normals,c_index,tx,ty,kx,ky,nx,ny,ncoef,n)


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
!     c_index - Integer, matrix size nx by ny: Index of coefficients for hte 
!               spline surface
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
  double precision, intent(in)     :: coef(ncoef,3)
  integer   , intent(in)           :: l_index(n,2)
  integer   , intent(in)           :: g_index(n)
  double precision, intent(out)    :: normals(n,3)
  integer   , intent(in)           :: c_index(nx,ny)
  integer,    intent(in)           :: kx,ky
  double precision, intent(in)     :: tx(kx+nx)
  double precision, intent(in)     :: ty(ky+ny)

  ! Local
  double precision                 :: s_coef(nx,ny,3)
  double precision                 :: dx(3)
  double precision                 :: dy(3)
  double precision                 :: normal(3)
  double precision                 :: WORK(3*MAX(KX,KY)+KY)

  double precision b2val

!   print *,'In getctlnormals'
!   print *,'n:',n
!   print *,'ncoef:',ncoef
!   print *,'nx,ny:',nx,ny
  ! Copy to the coefficient matrix

  do i=1,nx
     do j =1,ny
        s_coef(i,j,:) = coef(c_index(i,j)+1,:)
     end do
  end do


  do icoef = 1,n
     
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
     normals(icoef,:) = -normal

  end do
end subroutine getctlnormals
