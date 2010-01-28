subroutine eval_surface(u,v,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surface evaluates the n-dimensional b-spline surface
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, u coordinate
  !     v       - Real, v coordinate
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     coef    - Real,Array of B-spline coefficients  Size (nctlu,nctlv,ndim)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     ndim    - Integer, Spatial Dimension
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim
  
  
  ! Input
  integer         , intent(in)          :: ku,kv,nctlu,nctlv,ndim
  double precision, intent(in)          :: u,v
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
  double precision, intent(in)          :: coef(nctlu,nctlv,ndim)

  ! Output
  double precision, intent(out)         :: val(ndim)

  ! Working
  integer                               :: idim
  double precision                      :: work(3*max(ku,kv) + kv)

  double precision b2val

  do idim=1,ndim
     val(idim) = b2val(u,v,0,0,tu,tv,nctlu,nctlv,ku,kv,coef(:,:,idim),work)
  end do
  
end subroutine eval_surface

subroutine eval_surface_V(u,v,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,n,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surface_V is a vector version of eval_surface
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, u coordinate, length(n)
  !     v       - Real, v coordinate, length(n)
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     coef    - Real,Array of B-spline coefficients  Size (nctlu,nctlv,ndim)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     ndim    - Integer, Spatial Dimension
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size n,ndim
  
  
  ! Input
  integer         , intent(in)          :: ku,kv,nctlu,nctlv,ndim,n
  double precision, intent(in)          :: u(n),v(n)
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
  double precision, intent(in)          :: coef(nctlu,nctlv,ndim)

  ! Output
  double precision, intent(out)         :: val(n,ndim)

  ! Working
  integer                               :: idim,i
  double precision                      :: work(3*max(ku,kv) + kv)

  double precision b2val

  do i=1,n
     do idim=1,ndim
        val(i,idim) = b2val(u(i),v(i),0,0,tu,tv,nctlu,nctlv,ku,kv,coef(:,:,idim),work)
     end do
  end do
end subroutine eval_surface_V

subroutine eval_surface_M(u,v,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,n,m,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surface_M is a matrix version of eval_surface
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, u coordinate, size(n,m)
  !     v       - Real, v coordinate, size(n,m)
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     coef    - Real,Array of B-spline coefficients  Size (nctlu,nctlv,ndim)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     ndim    - Integer, Spatial Dimension
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size (n,m,ndim)
  
  
  ! Input
  integer         , intent(in)          :: ku,kv,nctlu,nctlv,ndim,n,m
  double precision, intent(in)          :: u(n,m),v(n,m)
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
  double precision, intent(in)          :: coef(nctlu,nctlv,ndim)

  ! Output
  double precision, intent(out)         :: val(n,m,ndim)

  ! Working
  integer                               :: idim,i,j
  double precision                      :: work(3*max(ku,kv) + kv)

  double precision b2val

  do i=1,n
     do j=1,m
        do idim=1,ndim
           val(i,j,idim) = b2val(u(i,j),v(i,j),0,0,tu,tv,nctlu,nctlv,ku,kv,coef(:,:,idim),work)
        end do
     end do
  end do
end subroutine eval_surface_M
