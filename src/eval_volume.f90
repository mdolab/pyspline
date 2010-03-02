subroutine eval_volume(u,v,w,tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surface evaluates the n-dimensional b-spline surface
  !              (SCLAR VERSION)
  !     Description of Arguments
  !     Input
  !     u       - Real, u coordinate
  !     v       - Real, v coordinate
  !     w       - Real, w coordinate
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     tw      - Real,Knot vector in w. Length nctlv+kw
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     kw      - Integer, order of B-spline in w
  !     coef    - Real,Array of B-spline coefficients  Size (nctlu,nctlv,nctlw,ndim)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     nctlw   - Integer,Number of control points in w
  !     ndim    - Integer, Spatial Dimension
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim
  
  
  ! Input
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,ndim
  double precision, intent(in)          :: u,v,w
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  double precision, intent(in)          :: coef(nctlu,nctlv,nctlw,ndim)

  ! Output
  double precision, intent(out)         :: val(ndim)

  ! Working
  integer                               :: idim
  double precision                      :: work(kv*kw+3*max(ku,kv,kw)+kw)

  double precision b3val

  do idim=1,ndim
     val(idim) = b3val(u,v,w,0,0,0,tu,tv,tw,nctlu,nctlv,nctlw,ku,kv,kw,coef(:,:,:,idim),work)
  end do
  
end subroutine eval_volume

subroutine eval_volume_V(u,v,w,tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,n,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surface evaluates the n-dimensional b-spline volume
  !              (VECTOR VERSION)
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, u coordinate, length(n)
  !     v       - Real, v coordinate, length(n)
  !     w       - Real, w coordinate, length(n)
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     tw      - Real,Knot vector in w. Length nctlv+kw
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     kw      - Integer, order of B-spline in w
  !     coef    - Real,Array of B-spline coefficients  Size (nctlu,nctlv,nctlw,ndim)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     nctlw   - Integer,Number of control points in w
  !     ndim    - Integer, Spatial Dimension
  !
  !     Ouput 
  !     val     - Real, Evaluated poinst, size (n,ndim)
  
  
  ! Input
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,ndim,n
  double precision, intent(in)          :: u(n),v(n),w(n)
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  double precision, intent(in)          :: coef(nctlu,nctlv,nctlw,ndim)

  ! Output
  double precision, intent(out)         :: val(n,ndim)

  ! Working
  integer                               :: idim,i,j
  double precision                      :: work(kv*kw+3*max(ku,kv,kw)+kw)

  double precision b3val
  do i=1,n
     do idim=1,ndim
        val(i,idim) = b3val(u(i),v(i),w(i),0,0,0,tu,tv,tw,nctlu,nctlv,nctlw,ku,kv,kw,coef(:,:,:,idim),work)
     end do
  end do
  
end subroutine eval_volume_V


subroutine eval_volume_M(u,v,w,tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,n,m,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surface_M evaluates the n-dimensional b-spline volume
  !              (MATRIX VERSION)
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, u coordinate, size(n,m)
  !     v       - Real, v coordinate, size(n,m)
  !     w       - Real, w coordinate, size(n,m)
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     tw      - Real,Knot vector in w. Length nctlv+kw
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     kw      - Integer, order of B-spline in w
  !     coef    - Real,Array of B-spline coefficients  Size (nctlu,nctlv,nctlw,ndim)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     nctlw   - Integer,Number of control points in w
  !     ndim    - Integer, Spatial Dimension
  !
  !     Ouput 
  !     val     - Real, Evaluated points, size (n,m,ndim)
  
  
  ! Input
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,ndim,n,m
  double precision, intent(in)          :: u(n,m),v(n,m),w(n,m)
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  double precision, intent(in)          :: coef(nctlu,nctlv,nctlw,ndim)

  ! Output
  double precision, intent(out)         :: val(n,m,ndim)

  ! Working
  integer                               :: idim,i,j
  double precision                      :: work(kv*kw+3*max(ku,kv,kw)+kw)

  double precision b3val
  do i=1,n
     do j=1,m
        do idim=1,ndim
           val(i,j,idim) = b3val(u(i,j),v(i,j),w(i,j),0,0,0,tu,tv,tw,nctlu,nctlv,nctlw,ku,kv,kw,coef(:,:,:,idim),work)
        end do
     end do
  end do
  
end subroutine eval_volume_M

subroutine eval_volume_T(u,v,w,tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,n,m,l,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surface_M evaluates the n-dimensional b-spline volume
  !              (TENSOR VERSION)
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, u coordinate, size(n,m,l)
  !     v       - Real, v coordinate, size(n,m,l)
  !     w       - Real, w coordinate, size(n,m,l)
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     tw      - Real,Knot vector in w. Length nctlv+kw
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     kw      - Integer, order of B-spline in w
  !     coef    - Real,Array of B-spline coefficients  Size (nctlu,nctlv,nctlw,ndim)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     nctlw   - Integer,Number of control points in w
  !     ndim    - Integer, Spatial Dimension
  !
  !     Ouput 
  !     val     - Real, Evaluated points, size (n,m,l,ndim)
  
  
  ! Input
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,ndim,n,m,l
  double precision, intent(in)          :: u(n,m,l),v(n,m,l),w(n,m,l)
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  double precision, intent(in)          :: coef(nctlu,nctlv,nctlw,ndim)

  ! Output
  double precision, intent(out)         :: val(n,m,l,ndim)

  ! Working
  integer                               :: idim,i,j,k
  double precision                      :: work(kv*kw+3*max(ku,kv,kw)+kw)

  double precision b3val
  do i=1,n
     do j=1,m
        do k=1,l
           do idim=1,ndim
              val(i,j,k,idim) = b3val(u(i,j,k),v(i,j,k),w(i,j,k),0,0,0,tu,tv,tw,nctlu,nctlv,nctlw,ku,kv,kw,coef(:,:,:,idim),work)
           end do
        end do
     end do
  end do
  
end subroutine eval_volume_T

subroutine eval_volume_deriv(u,v,w,tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surface evaluates the n-dimensional b-spline volume derivative
  !              (SCLAR VERSION)
  !     Description of Arguments
  !     Input
  !     u       - Real, u coordinate
  !     v       - Real, v coordinate
  !     w       - Real, w coordinate
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     tw      - Real,Knot vector in w. Length nctlv+kw
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     kw      - Integer, order of B-spline in w
  !     coef    - Real,Array of B-spline coefficients  Size (nctlu,nctlv,nctlw,ndim)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     nctlw   - Integer,Number of control points in w
  !     ndim    - Integer, Spatial Dimension
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim
  
  
  ! Input
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,ndim
  double precision, intent(in)          :: u,v,w
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  double precision, intent(in)          :: coef(nctlu,nctlv,nctlw,ndim)

  ! Output
  double precision, intent(out)         :: val(3,ndim)

  ! Working
  integer                               :: idim
  double precision                      :: work(kv*kw+3*max(ku,kv,kw)+kw)

  double precision b3val

  do idim=1,ndim
     val(1,idim) = b3val(u,v,w,1,0,0,tu,tv,tw,nctlu,nctlv,nctlw,ku,kv,kw,coef(:,:,:,idim),work)
     val(2,idim) = b3val(u,v,w,0,1,0,tu,tv,tw,nctlu,nctlv,nctlw,ku,kv,kw,coef(:,:,:,idim),work)
     val(3,idim) = b3val(u,v,w,0,0,1,tu,tv,tw,nctlu,nctlv,nctlw,ku,kv,kw,coef(:,:,:,idim),work)
  end do
  
end subroutine eval_volume_deriv


subroutine eval_volume_deriv2(u,v,w,tu,tv,tw,ku,kv,kw,coef,nctlu,nctlv,nctlw,ndim,val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surface evaluates the n-dimensional b-spline volume second derivative
  !              (SCLAR VERSION)
  !     Description of Arguments
  !     Input
  !     u       - Real, u coordinate
  !     v       - Real, v coordinate
  !     w       - Real, w coordinate
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     tw      - Real,Knot vector in w. Length nctlv+kw
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     kw      - Integer, order of B-spline in w
  !     coef    - Real,Array of B-spline coefficients  Size (nctlu,nctlv,nctlw,ndim)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     nctlw   - Integer,Number of control points in w
  !     ndim    - Integer, Spatial Dimension
  !
  !     Ouput 
  !     val     - Real, Evaluated point, size ndim
  
  
  ! Input
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,ndim
  double precision, intent(in)          :: u,v,w
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  double precision, intent(in)          :: coef(nctlu,nctlv,nctlw,ndim)

  ! Output
  double precision, intent(out)         :: val(3,3,ndim)

  ! Working
  integer                               :: idim
  double precision                      :: work(kv*kw+3*max(ku,kv,kw)+kw)

  double precision b3val

  do idim=1,ndim

     ! Row 1
     if (ku>=3) then
        val(1,1,idim) = b3val(u,v,w,2,0,0,tu,tv,tw,nctlu,nctlv,nctlw,ku,kv,kw,coef(:,:,:,idim),work)
     else
        val(1,1,idim) = 0.0
     end if
     
     val(1,2,idim) = b3val(u,v,w,1,1,0,tu,tv,tw,nctlu,nctlv,nctlw,ku,kv,kw,coef(:,:,:,idim),work)
     val(1,3,idim) = b3val(u,v,w,1,0,1,tu,tv,tw,nctlu,nctlv,nctlw,ku,kv,kw,coef(:,:,:,idim),work)

     ! Row 2
     val(2,1,idim) = val(1,2,dim)
     if (kv>=3) then
        val(2,2,idim) = b3val(u,v,w,0,2,0,tu,tv,tw,nctlu,nctlv,nctlw,ku,kv,kw,coef(:,:,:,idim),work)
     else
        val(2,2,idim) = 0.0
     end if
     val(2,3,idim) = b3val(u,v,w,0,1,1,tu,tv,tw,nctlu,nctlv,nctlw,ku,kv,kw,coef(:,:,:,idim),work)

     ! Row 3
     val(3,1,idim) = val(1,3,idim)
     val(3,2,idim) = val(2,3,idim)
     if (kw>=3) then
        val(3,3,idim) = b3val(u,v,w,0,0,2,tu,tv,tw,nctlu,nctlv,nctlw,ku,kv,kw,coef(:,:,:,idim),work)
     else
        val(3,3,idim) = 0.0
     end if
  end do
  
end subroutine eval_volume_deriv2
