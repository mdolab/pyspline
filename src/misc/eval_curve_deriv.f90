subroutine eval_curve_deriv(s,t,k,coef,nctl,ndim,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surf_deriv is a scalar version of the 
  !              NURBS derivative evaluation function
  !
  !     Description of Arguments
  !     Input
  !     s       - Real, s coordinate
  !     t       - Real,Knot vector. Length nctl+k
  !     k       - Integer,order of NURBS
  !     coef    - Real,Array of NURBS coefficients and weights. Size (nctl,ndim+1)
  !     nctl   - Integer,Number of control points
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim
  
  ! Input
  integer         , intent(in)          :: k,nctl,ndim
  double precision, intent(in)          :: s
  double precision, intent(in)          :: t(nctl+k)
  double precision, intent(in)          :: coef(nctl,ndim+1)

  ! Output
  double precision, intent(out)         :: val(ndim)

  ! Working
  integer                               :: idim,inbv
  double precision                      :: work(3*k),deriv(ndim)
  double precision                      :: weight

  !double precision bvalu

  inbv = 1

  ! Evaluate the weight first
  weight = bvalu(t,coef(:,ndim+1),nctl,k,0,s,inbv,work)
  weight_deriv = bvalu(t,coef(:,ndim+1),nctl,k,1,s,inbv,work)
  do idim=1,ndim
     val(idim) = bvalu(t,coef(:,idim)*coef(:,ndim+1),nctl,k,0,s,inbv,work)
     deriv(idim) = bvalu(t,coef(:,idim)*coef(:,ndim+1),nctl,k,1,s,inbv,work)
  end do
  
  val = (deriv-weight_deriv*val)/weight
  
end subroutine eval_curve_deriv
