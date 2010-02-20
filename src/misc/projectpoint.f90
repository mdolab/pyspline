subroutine projectpoint(coef,kx,ky,nx,ny,tx,ty,x0,u0,v0,Niter,tol,D,converged)


!*** DESCRIPTION
!
!     Written by Gaetan Kenway
! 
!     The function takes the spline defination of a bivariate spline
!     surface as defined by coef,kx,ky,nx,ny,tx,ty and a point x0 and
!     determines the u,v positions that minimizes the distance between
!     the point and the surface. 

!     Description of Arguments:
!     Input:
!     coef    - Real, Matrix of coefficients, size nx by ny by 3 
!               (This only works with 3 spatial degrees of freedom)
!     kx      - Integer, The spline order in x
!     ky      - Integer, The spline order in y
!     nx      - Integer: Number of contol points in x
!     ny      - Integer: Number of control points in y
!     tx      - Real, vector of knots in x, size nx + kx
!     ty      - Real, vector of knots in y, size ny + ky
!     x0      - Real, vector or length 3: Point to project onto the surface
!     u0      - Real: Initial guess for u
!     v0      - Real: Initial guess for v
!     Niter   - Integer: Maximum number of iterations
!     tol     - Real: Tolerance for newton iteration

!     Output:
!     u0      - Real: Surface Position of point projection
!     v0      - Real: Surface Position of point projection
!     D       - Real: Distance between point and point projection
!     converged - Integer: 1 if converged -1 if not converged
!
  implicit none

  integer i,j,idim
  double precision, intent(in)     :: coef(nx,ny,3)
  integer         , intent(in)     :: kx,ky
  integer         , intent(in)     :: nx,ny
  double precision, intent(in)     :: tx(kx+nx)
  double precision, intent(in)     :: ty(ky+ny)
  double precision, intent(in)     :: X0(3)
  double precision, intent(inout)  :: u0,v0
  integer         , intent(in)     :: Niter
  double precision, intent(in)     :: tol

  double precision, intent(out)    :: D(3)
  integer         , intent(out)    :: converged

  double precision                 :: dx(3),dy(3),n(3)
  double precision                 :: u_length,v_length,updateu,updatev
  double precision                 :: D2(3)
  double precision                 :: D_norm,D2_norm
  double precision                 :: WORK(3*MAX(KX,KY)+KY+KX)
  double precision                 :: distance,base_point(3),norm,cur_point(3)
  double precision b2val

  converged = -1

!   print *,'Welcome to projectpoint!'
!   print *,'u0,v0:',u0,v0
!   print *,'nx,ny,kx,ky:',nx,ny,kx,ky
 
  do i =1,Niter

     do idim =1,3
        dx(idim) = b2val(u0,v0,1,0,tx,ty,nx,ny,kx,ky,coef(:,:,idim),work)
        dy(idim) = b2val(u0,v0,0,1,tx,ty,nx,ny,kx,ky,coef(:,:,idim),work)
     end do

     do idim =1,3
        cur_point(idim) = b2val(u0,v0,0,0,tx,ty,nx,ny,kx,ky,coef(:,:,idim),work)
     end do
     
      D = x0-cur_point

     ! We Need to project D onto dx and dy

     updateu = dot_product(D,dx)/(dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3))
     updatev = dot_product(D,dy)/(dy(1)*dy(1) + dy(2)*dy(2) + dy(3)*dy(3))

     if (u0 + updateu > tx(nx+Kx)) then
        u0 = tx(nx + kx)
        updateu = 0
     else if (u0 + updateu < tx(1)) then
        u0 = tx(1)
        updateu = 0
     else
        u0 = u0 + updateu
     end if

     if (v0 + updatev > ty(ny+ky)) then
        v0 = ty(ny + ky)
        updatev = 0
     else if (v0 + updatev < ty(1)) then
        v0 = ty(1)
        updatev = 0
     else
        v0 = v0 + updatev
     end if

     if (abs(updateu)<tol .and. abs(updatev)<tol) then
        converged = 1
        exit
     end if
        

  end do
  
end subroutine projectpoint
