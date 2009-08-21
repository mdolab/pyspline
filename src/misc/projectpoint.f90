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

  double precision                 :: dx(3),dy(3),n,norm
  double precision                 :: u_length,v_length,updateu,updatev
  double precision                 :: D2(3)
  double precision                 :: D_norm,D2_norm
  double precision                 :: WORK(3*MAX(KX,KY)+KY+KX)

  double precision b2val

  converged = -1

  !print *,'Welcome to projectpoint!'
 
  do i =1,Niter

     do idim =1,3
        dx(idim) = b2val(x0,y0,1,0,tx,ty,nx,ny,kx,ky,coef(:,:,idim),work)
        dy(idim) = b2val(x0,y0,0,1,tx,ty,nx,ny,kx,ky,coef(:,:,idim),work)
     end do

     ! Now get the cross product 

     n(1) = dx(2)*dy(3)-dx(3)*dy(2)
     n(2) = dx(3)*dy(1)-dx(1)*dy(3)
     n(3) = dx(1)*dy(2)-dx(2)*dy(1)

     ! Normalizet the cross product
     norm = sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))
     if (norm .ne. 0) then
        n = n /norm
     else
        print *,'The Normal Vector to the surface is ill-defined at'
        print *,'x = ',x0
        print *,'y = ',y0
        stop
     end if



     do idim =1,3
        D(idim) = x0(idim) - b2val(u0,v0,0,0,tx,ty,nx,ny,kx,ky,coef(:,:,idim),work)
     end d
     


     u_length = sqrt(Ydotu(1)*Ydotu(1) + Ydotu(2)*Ydotu(2) + Ydotu(3)*Ydotu(3))
     v_length = sqrt(Ydotv(1)*Ydotv(1) + Ydotv(2)*Ydotv(2) + Ydotv(3)*Ydotv(3))

     updateu = 0
     updatev = 0
     do idim = 1,3
        if (Ydotu(idim) .ne. 0 .and. u_length .ne. 0) then
           updateu = updateu + D(idim)*Ydotu(idim)/u_length
        end if
        if (Ydotv(idim) .ne. 0 .and. v_length .ne. 0) then
           updatev = updatev + D(idim)*Ydotv(idim)/v_length
        end if
     end do
           
     if (u0+updateu > tx(nx+kx)) then
        updateu = tx(nx+kx) - u0
     else if (u0 + updateu < tx(1)) then
        updateu = tx(1) - u0
     end if
           
     if (v0+updatev > ty(ny+ky)) then
        updatev = ty(ny+ky) -v0
     else if (v0 + updatev < ty(1)) then
        updatev = ty(1) - v0
     end if

     do idim =1,3
        D2(idim) = x0(idim) - b2val(u0+updateu,v0+updatev,0,0,tx,ty,nx,ny,kx,ky,coef(:,:,idim),work)
     end do

     D2_norm = sqrt(D2(1)*D2(1) + D2(2)*D2(2) + D2(3)*D2(3))
     D_norm  = sqrt(D (1)*D (1) + D (2)*D (2) + D (3)*D (3))

     if (abs(D2_norm)>abs(D_norm)) then
        !  Half Update
        updateu = updateu/4
        updatev = updatev/4
     end if

     if (abs(updateu)<tol .and. abs(updatev)<tol) then
        u0 = u0 + updateu
        v0 = v0 + updatev
        converged = 1

        do idim =1,3 ! Get the final Difference
           D(idim) = x0(idim) - b2val(u0,v0,0,0,tx,ty,nx,ny,kx,ky,coef(:,:,idim),work)
        end do
      else
        u0 = u0 + updateu
        v0 = v0 + updatev
     end if
  end do
  
end subroutine projectpoint
