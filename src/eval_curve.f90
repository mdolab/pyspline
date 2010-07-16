subroutine eval_curve(s,t,k,coef,nctl,ndim,n,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surf is a vector version of the B-spline evaluation function
  !
  !     Description of Arguments
 !     Input
  !     s       - Real, Vector of s coordinates, length n
  !     t       - Real,Knot vector. Length nctl+k
  !     k       - Integer,order of B-spline 
  !     coef    - Real,Array of b-sline coefficients  Size (ndim,nctl)
  !     nctl    - Integer,Number of control points
  !
  !     Ouput 
  !     val     - Real, Evaluated points, size ndim by n
  implicit none
  ! Input
  integer         , intent(in)          :: k,nctl,ndim,n
  double precision, intent(in)          :: s(n)
  double precision, intent(in)          :: t(nctl+k)
  double precision, intent(in)          :: coef(ndim,nctl)

  ! Output
  double precision, intent(out)         :: val(ndim,n)

  ! Working
  integer                               :: i,l,idim,istart
  integer                               :: ILEFT,IWORK,ILO,mflag
  double precision                      :: basisu(k)

 val(:,:) = 0.0
 ILO = 1
  do i=1,n
     call INTRV(T,NCTL+K,s(i),ILO,ILEFT,MFLAG)
     if (mflag == 1) then
        ileft = ileft-k
     end if
     call basis(t,nctl,k,s(i),ileft,basisu)
     istart = ileft-k
     do l=1,k
        do idim=1,ndim
           val(idim,i) = val(idim,i) + basisu(l)*coef(idim,istart+l)
        end do
     end do
  end do
end subroutine eval_curve

subroutine eval_curve_deriv(s,t,k,coef,nctl,ndim,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surf_deriv is a scalar version of the 
  !              B-spline derivative evaluation function
  !
  !     Description of Arguments
  !     Input
  !     s       - Real, s coordinate
  !     t       - Real,Knot vector. Length nctl+k
  !     k       - Integer,order of B-spline
  !     coef    - Real,Array of B-spline coefficients. Size (ndim,nctl)
  !     nctl   - Integer,Number of control points
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim
  
  implicit none
  ! Input
  integer         , intent(in)          :: k,nctl,ndim
  double precision, intent(in)          :: s
  double precision, intent(in)          :: t(nctl+k)
  double precision, intent(in)          :: coef(ndim,nctl)

  ! Output
  double precision, intent(out)         :: val(ndim)

  ! Working
  integer                               :: idim,inbv
  double precision                      :: work(3*k)

  ! Functions
  double precision                      :: bvalu

  inbv = 1
  
  do idim=1,ndim
     val(idim) = bvalu(t,coef(idim,:),nctl,k,1,s,inbv,work)
  end do
    
end subroutine eval_curve_deriv

subroutine eval_curve_deriv2(s,t,k,coef,nctl,ndim,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_curve_deriv2 evals the 2nd derivative
  !              B-spline derivative evaluation function
  !
  !     Description of Arguments
  !     Input
  !     s       - Real, s coordinate
  !     t       - Real,Knot vector. Length nctl+k
  !     k       - Integer,order of B-spline
  !     coef    - Real,Array of B-spline coefficients Size (ndim,nctl)
  !     nctl    - Integer,Number of control points
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim
  
  implicit none
  ! Input
  integer         , intent(in)          :: k,nctl,ndim
  double precision, intent(in)          :: s
  double precision, intent(in)          :: t(nctl+k)
  double precision, intent(in)          :: coef(ndim,nctl)

  ! Output
  double precision, intent(out)         :: val(ndim)

  ! Working
  integer                               :: idim,inbv
  double precision                      :: work(3*k)

  ! Functions
  double precision                      :: bvalu

  inbv = 1
  
  if (k == 2) then
     val(:) = 0.0
  else
     do idim=1,ndim
        val(idim) = bvalu(t,coef(idim,:),nctl,k,2,s,inbv,work)
     end do
  end if

end subroutine eval_curve_deriv2


subroutine eval_curve_c(s,t,k,coef,nctl,ndim,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_curve_complex is a special version for evaluation
  !     of data when control points are used in a complex-step analysis
  !     The actual derivative is computed analytically, but the output looks
  !     like the CS went all the way through
  !
  !     Description of Arguments
  !     Input
  !     s       - Real, s coordinate
  !     t       - Real,Knot vector. Length nctl+k
  !     k       - Integer,order of B-spline
  !     coef    - Real,Array of B-spline coefficients  Size (ndim,nctl)
  !     nctl   - Integer,Number of control points
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim
  implicit none
  ! Input
  integer         , intent(in)          :: k,nctl,ndim
  double precision, intent(in)          :: s
  double precision, intent(in)          :: t(nctl+k)
  complex*16      , intent(in)          :: coef(ndim,nctl)

  ! Output
  complex*16      , intent(out)         :: val(ndim)

  ! Working
  integer                               :: idim,ilow,ileft,mflag,iwork,i,ilo
  double precision                      :: work(3*k),temp_val,temp_deriv
  double precision                      :: basisu(k)

  ! Functions
  double precision bvalu

  
  ilo = 1

  do idim=1,ndim
     call intrv(t,nctl+k,s,ilo,ileft,mflag)
     if (mflag == 0) then
        call basis(t,nctl,k,s,ileft,basisu)
     else if (mflag == 1) then
        ileft = nctl
        basisu(:) = 0.0
        basisu(k) = 1.0
     end if

     ! Get the value...we can use the ileft computed above
     temp_val = bvalu(t,real(coef(idim,:)),nctl,k,0,s,ileft,work)
     ! Now get the derivative
     temp_deriv = 0.0
     do i=1,k
        temp_deriv = temp_deriv + basisu(i)*aimag(coef(idim,ileft-k+i))
     end do

     val(idim) = cmplx(temp_val,temp_deriv)

  end do
  
end subroutine eval_curve_c
