subroutine eval_surface(u, v, tu, tv, ku, kv, coef, nctlu, nctlv, ndim, &
     n, m, val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surface evaluates (possibly) many points on the surface
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, u coordinate, size(m,n)
  !     v       - Real, v coordinate, size(m,n)
  !     tu      - Real, Knot vector in u. size(nctlu+ku)
  !     tv      - Real, Knot vector in v. size(nctlv+kv)
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     coef    - Real, Array of B-spline coefficients  Size (ndim,nctlv,nctlu)
  !     nctlu   - Integer, Number of control points in u
  !     nctlv   - Integer, Number of control points in v
  !     ndim    - Integer, Spatial Dimension
  !
  !     Ouput 
  !     val     - Real, Evaluated point(s), size (ndim,m,n)

  implicit none

  ! Input
  integer         , intent(in)    :: ku,kv,nctlu,nctlv,ndim,n,m
  double precision, intent(in)    :: u(m,n),v(m,n)
  double precision, intent(in)    :: tu(nctlu+ku),tv(nctlv+kv)
  double precision, intent(in)    :: coef(ndim,nctlv,nctlu)

  ! Output
  double precision, intent(out)   :: val(ndim,m,n)

  ! Working
  integer                         :: idim,istartu,istartv,i,j,ii,jj
  integer                         :: ileftu,ilou,mflagu
  integer                         :: ileftv,ilov,mflagv
  double precision                :: basisu(ku),basisv(kv)

  ilou = 1
  ilov = 1
  val(:,:,:) = 0.0
  do ii=1,n
     do jj = 1,m
        ! U
        call INTRV(tu,nctlu+ku,u(jj,ii),ilou,ileftu,mflagu)
        if (mflagu == 1) then
           ileftu = ileftu-ku
        end if
        call basis(tu,nctlu,ku,u(jj,ii),ileftu,basisu)
        istartu = ileftu-ku

        ! V
        call INTRV(tv,nctlv+kv,v(jj,ii),ilov,ileftv,mflagv)
        if (mflagv == 1) then
           ileftv = ileftv-kv
        end if
        call basis(tv,nctlv,kv,v(jj,ii),ileftv,basisv)
        istartv = ileftv-kv

        do i=1,ku
           do j=1,kv
              do idim=1,ndim
                 val(idim,jj,ii) = val(idim,jj,ii) + basisu(i)*basisv(j)*coef(idim,istartv+j,istartu+i)
              end do
           end do
        end do
     end do
  end do

end subroutine eval_surface

subroutine eval_surface_deriv(u, v, tu, tv, ku, kv, coef, nctlu, nctlv, &
     ndim, val)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surface_deriv evaluates the derivative of a
  !     point on the surface
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, u coordinate
  !     v       - Real, v coordinate
  !     tu      - Real, Knot vector in u. size(nctlu+ku)
  !     tv      - Real, Knot vector in v. size(nctlv+kv)
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     coef    - Real, Array of B-spline coefficients  Size (ndim,nctlv,nctlu)
  !     nctlu   - Integer, Number of control points in u
  !     nctlv   - Integer, Number of control points in v
  !     ndim    - Integer, Spatial Dimension
  !
  !     Ouput 
  !     val     - Real, Evaluated derivatives, size (ndim,2)

  ! Input
  integer         , intent(in)          :: ku,kv,nctlu,nctlv,ndim
  double precision, intent(in)          :: u,v
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
  double precision, intent(in)          :: coef(ndim,nctlv,nctlu)

  ! Output
  double precision, intent(out)         :: val(ndim,2)

  ! Working
  integer                               :: idim
  double precision                      :: work(3*max(ku,kv) + kv)
  double precision b2val

  do idim=1,ndim
     val(idim,1) = b2val(v,u,0,1,tv,tu,nctlv,nctlu,kv,ku,coef(idim,:,:),work)
     val(idim,2) = b2val(v,u,1,0,tv,tu,nctlv,nctlu,kv,ku,coef(idim,:,:),work)
  end do

end subroutine eval_surface_deriv

subroutine eval_surface_deriv2(u,v,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)


  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract eval_surface_deriv2 evaluates the second derivative
  !     of a point on the surface
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, u coordinate
  !     v       - Real, v coordinate
  !     tu      - Real, Knot vector in u. size(nctlu+ku)
  !     tv      - Real, Knot vector in v. size(nctlv+kv)
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     coef    - Real, Array of B-spline coefficients  Size (ndim,nctlv,nctlu)
  !     nctlu   - Integer, Number of control points in u
  !     nctlv   - Integer, Number of control points in v
  !     ndim    - Integer, Spatial Dimension
  !
  !     Ouput 
  !     val     - Real, Evaluated second derivatives, size (ndim,2,2)

  ! Input
  integer         , intent(in)          :: ku,kv,nctlu,nctlv,ndim
  double precision, intent(in)          :: u,v
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
  double precision, intent(in)          :: coef(ndim,nctlv,nctlv)

  ! Output
  double precision, intent(out)         :: val(ndim,2,2)

  ! Working
  integer                               :: idim
  double precision                      :: work(3*max(ku,kv) + kv)

  double precision b2val

  do idim=1,ndim
     if (ku>=3) then
        val(idim,1,1) = b2val(v,u,0,2,tv,tu,nctlv,nctlu,kv,ku,coef(idim,:,:),work)
     else
        val(idim,1,1) = 0.0
     end if

     val(idim,1,2) = b2val(v,u,1,1,tv,tu,nctlv,nctlu,kv,ku,coef(idim,:,:),work)
     val(idim,2,1) = val(idim,1,2)

     if (kv>=3) then
        val(idim,2,2) = b2val(v,u,2,0,tv,tu,nctlv,nctlu,kv,ku,coef(idim,:,:),work)
     else
        val(idim,2,2) = 0.0
     end if
  end do

end subroutine eval_surface_deriv2
