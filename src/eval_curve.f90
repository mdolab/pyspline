subroutine eval_curve(s,t,k,coef,nctl,ndim,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_curve is a scalar version of the B-spline evaluation function
  !
  !     Description of Arguments
  !     Input
  !     s       - Real, s coordinate
  !     t       - Real,Knot vector. Length nctl+k
  !     k       - Integer,order of B-spline
  !     coef    - Real,Array of B-spline coefficients  Size (nctl,ndim)
  !     nctl   - Integer,Number of control points
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim
  implicit none
  ! Input
  integer         , intent(in)          :: k,nctl,ndim
  double precision, intent(in)          :: s
  double precision, intent(in)          :: t(nctl+k)
  double precision, intent(in)          :: coef(nctl,ndim)

  ! Output
  double precision, intent(out)         :: val(ndim)

  ! Working
  integer                               :: idim,inbv
  double precision                      :: work(3*k)

  ! Functions
  double precision bvalu

  inbv = 1

  do idim=1,ndim
     val(idim) = bvalu(t,coef(:,idim),nctl,k,0,s,inbv,work)
  end do
  
end subroutine eval_curve

subroutine eval_curve_V(s,t,k,coef,nctl,ndim,n,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surf_V is a vector version of the B-spline evaluation function
  !
  !     Description of Arguments
 !     Input
  !     s       - Real, Vector of s coordinates, length n
  !     t       - Real,Knot vector. Length nctl+k
  !     k       - Integer,order of B-spline 
  !     coef    - Real,Array of b-sline coefficients  Size (nctl,ndim)
  !     nctl    - Integer,Number of control points
  !
  !     Ouput 
  !     val     - Real, Evaluated points, size n by ndim
  implicit none
  ! Input
  integer         , intent(in)          :: k,nctl,ndim,n
  double precision, intent(in)          :: s(n)
  double precision, intent(in)          :: t(nctl+k)
  double precision, intent(in)          :: coef(nctl,ndim)

  ! Output
  double precision, intent(out)         :: val(n,ndim)

  ! Working
  integer                               :: i,idim,inbv
  double precision                      :: work(3*k)

  ! Functions
  double precision                      :: bvalu

  inbv = 1

  do i=1,n
     do idim=1,ndim
        val(i,idim) = bvalu(t,coef(:,idim),nctl,k,0,s(i),inbv,work)
     end do
  end do

end subroutine eval_curve_V
