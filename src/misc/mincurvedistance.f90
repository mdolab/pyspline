subroutine mincurvedistance(t1,k1,coef1,t2,k2,coef2,n1,n2,s,t,Niter,tol,D,converged)

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
!     t1      - Knot Vector for Curve 1
!     k1      - Order for Curve 1
!     coef1   - Coefficients for Curve 1
!     t2      - Knot Vector for Curve 2
!     k2      - Order for Curve 2
!     coef2   - Coefficients for Curve 2
!     n1      - Number of coefficients on curve 1
!     n2      - Number of coefficients on curve 2
!     Niter   - Integer: Maximum number of iterations
!     tol     - Real: Tolerance for newton iteration
!     s       - Real: Initial guess for parameter on Curve 1
!     t       - Real: Initial guess for parameter on Curve 2

!     Output:
!     s       - Real: parameter on Curve 1
!     t       - Real: parameter on Curve 2
!     D       - Real: Distance between curve 
!     converged - Integer: 1 if converged -1 if not converged
!
  implicit none

  integer i,idim
  integer         , intent(in)     :: n1,n2,k1,k2
  double precision, intent(in)     :: t1(n1+k1),t2(n2+k2)
  double precision, intent(in)     :: coef1(n1,3),coef2(n2,3)
  double precision, intent(inout)  :: s,t 
  integer         , intent(in)     :: Niter
  double precision, intent(in)     :: tol

  double precision, intent(out)    :: D
  integer         , intent(out)    :: converged

  double precision                 :: WORK1(3*k1),WORK2(3*k2)
  integer                          :: INBV1,INBV2
  double precision                 :: D1(3),D2(3),D3
  double precision                 :: dDdu(2)
  double precision                 :: J(2,2),invJ(2,2)
  double precision                 ::   x1(3),  x2(3)
  double precision                 ::  dx1(3), dx2(3)
  double precision                 :: ddx1(3),ddx2(3)
  double precision                 :: inv_det
  double precision                 :: update(2)
  integer                          :: sub_iter
  

  double precision bvalu

  converged = -1

!   print *,'Welcome to minCurveDistance!'
!   print *,'k1,n1,k2,n2',k1,n1,k2,n2
!   print *,'t1,t2:',t1,t2
!   print *,'coef1:',coef1
  inbv1 = 1
  inbv2 = 1
  sub_iter = 1
  do i =1,Niter
     !print *,'iter,s,t:',i,s,t,sub_iter

     do idim =1,3
        x1(idim) = bvalu(t1,coef1(:,idim),n1,k1,0,s,INBV1,WORK1)
        x2(idim) = bvalu(t2,coef2(:,idim),n2,k2,0,t,INBV2,WORK2)
      
        dx1(idim) = bvalu(t1,coef1(:,idim),n1,k1,1,s,INBV1,WORK1)
        dx2(idim) = bvalu(t2,coef2(:,idim),n2,k2,1,t,INBV2,WORK2)
      
        if (k1 > 2) then
           ddx1(idim) = bvalu(t1,coef1(:,idim),n1,k1,2,s,INBV1,WORK1)
        else
           ddx1(idim) = 0
        end if
        if (k2 > 2) then
           ddx2(idim) = bvalu(t2,coef2(:,idim),n2,k2,2,t,INBV2,WORK2)
        else
           ddx2(idim) = 0
        end if
     end do

     D1 = x1-x2
     D = D1(1)*D1(1) + D1(2)*D1(2) + D1(3)*D1(3)
     ! D is the squared distance we are minimizing

     dDdu(1) =  dot_product(D1,dx1)
     dDdu(2) = -dot_product(D1,dx2)

     J(1,1) = dot_product(D1,ddx1) + dot_product(dx1,dx1)
     J(1,2) = -dot_product(dx2,dx1)
     J(2,1) = J(1,2)
     J(2,2) = dot_product(D1,ddx2) + dot_product(dx2,dx2)

     inv_det = 1/(J(1,1)*J(2,2)-J(1,2)*J(2,1))

     invJ(1,1) = inv_det * J(2,2)
     invJ(1,2) = -inv_det * J(2,1)
     invJ(2,1) = -inv_det * J(1,2)
     invJ(2,2) = inv_det * J(1,1)

     update = matmul(invJ,dDdu)

     sub_iter = 1
     do while(sub_iter < 10)
        !print *,'sub_iter:',sub_iter
        ! Check to see if we have a bad update
        do idim =1,3
           x1(idim) = bvalu(t1,coef1(:,idim),n1,k1,0,s-update(1),INBV1,WORK1)
           x2(idim) = bvalu(t2,coef2(:,idim),n2,k2,0,t-update(2),INBV2,WORK2)
        end do
        D2 = x1-x2
        D3 = D2(1)*D2(1) + D2(2)*D2(2) + D2(3)*D2(3)
        !print *,'D3:',D3
        if (D3 < D) then ! Accept the step and stop
           s = s - update(1)
           t = t - update(2)
           exit
        else if (D3 > D) then
           update = update/2
        end if

        if (sub_iter .eq. 9) then
           ! we are at end...take update anyway
           s = s - update(1)
           t = t - update(2)
        end if
        
        sub_iter = sub_iter + 1
     end do

     if (abs(update(1))<tol .and. abs(update(2))<tol) then
        converged = 1
        exit
     end if
  end do
  
end subroutine mincurvedistance
