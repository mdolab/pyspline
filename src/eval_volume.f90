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
  
  implicit none
  ! Input
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,ndim
  double precision, intent(in)          :: u,v,w
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  double precision, intent(in)          :: coef(nctlu,nctlv,nctlw,ndim)

  ! Output
  double precision, intent(out)         :: val(ndim)

 ! Working
  integer                               :: idim,istartu,istartv,istartw,i,j,k
  integer                               :: ileftu,iworku,ilou,mflagu
  integer                               :: ileftv,iworkv,ilov,mflagv
  integer                               :: ileftw,iworkw,ilow,mflagw
  double precision                      :: basisu(ku),basisv(kv),basisw(kw)
  double precision                      :: worku(2*ku),workv(2*kv),workw(2*kw)

  val(:) = 0.0
  ilou = 1
  ilov = 1
  ilow = 1
  ! U
  call INTRV(tu,nctlu+ku,u,ilou,ileftu,mflagu)
  if (mflagu == 1) then
     ileftu = ileftu-ku
  end if
  call BSPVN(tu,ku,ku,1,u,ileftu,basisu,worku,iworku)
  istartu = ileftu-ku
  
  ! V
  call INTRV(tv,nctlv+kv,v,ilov,ileftv,mflagv)
  if (mflagv == 1) then
     ileftv = ileftv-kv
  end if
  call BSPVN(tv,kv,kv,1,v,ileftv,basisv,workv,iworkv)
  istartv = ileftv-kv

  ! W
  call INTRV(tw,nctlw+kw,w,ilow,ileftw,mflagw)
  if (mflagw == 1) then
     ileftw = ileftw-kw
  end if
  call BSPVN(tw,kw,kw,1,w,ileftw,basisw,workw,iworkw)
  istartw = ileftw-kw

  do i=1,ku
     do j=1,kv
        do k=1,kw
           do idim=1,ndim
              val(idim) = val(idim) + basisu(i)*basisv(j)*basisw(k)*coef(istartu+i,istartv+j,istartw+k,idim)
           end do
        end do
     end do
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
  
  implicit none
  ! Input
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,ndim,n
  double precision, intent(in)          :: u(n),v(n),w(n)
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  double precision, intent(in)          :: coef(nctlu,nctlv,nctlw,ndim)

  ! Output
  double precision, intent(out)         :: val(n,ndim)

! Working
  integer                               :: idim,istartu,istartv,istartw,i,j,k,ii
  integer                               :: ileftu,iworku,ilou,mflagu
  integer                               :: ileftv,iworkv,ilov,mflagv
  integer                               :: ileftw,iworkw,ilow,mflagw
  double precision                      :: basisu(ku),basisv(kv),basisw(kw)
  double precision                      :: worku(2*ku),workv(2*kv),workw(2*kw)

  val(:,:) = 0.0
  ilou = 1
  ilov = 1
  ilow = 1
  
  do ii=1,n
     ! U
     call INTRV(tu,nctlu+ku,u(ii),ilou,ileftu,mflagu)
     if (mflagu == 1) then
        ileftu = ileftu-ku
     end if
     call BSPVN(tu,ku,ku,1,u(ii),ileftu,basisu,worku,iworku)
     istartu = ileftu-ku
     
     ! V
     call INTRV(tv,nctlv+kv,v(ii),ilov,ileftv,mflagv)
     if (mflagv == 1) then
        ileftv = ileftv-kv
     end if
     call BSPVN(tv,kv,kv,1,v(ii),ileftv,basisv,workv,iworkv)
     istartv = ileftv-kv
     
     ! W
     call INTRV(tw,nctlw+kw,w(ii),ilow,ileftw,mflagw)
     if (mflagw == 1) then
        ileftw = ileftw-kw
     end if
     call BSPVN(tw,kw,kw,1,w(ii),ileftw,basisw,workw,iworkw)
     istartw = ileftw-kw

     do i=1,ku
        do j=1,kv
           do k=1,kw
              do idim=1,ndim
                 val(ii,idim) = val(ii,idim) + basisu(i)*basisv(j)*basisw(k)*coef(istartu+i,istartv+j,istartw+k,idim)
              end do
           end do
        end do
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
  
  implicit none
  ! Input
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,ndim,n,m
  double precision, intent(in)          :: u(n,m),v(n,m),w(n,m)
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  double precision, intent(in)          :: coef(nctlu,nctlv,nctlw,ndim)

  ! Output
  double precision, intent(out)         :: val(n,m,ndim)

! Working
  integer                               :: idim,istartu,istartv,istartw,i,j,k,ii,jj
  integer                               :: ileftu,iworku,ilou,mflagu
  integer                               :: ileftv,iworkv,ilov,mflagv
  integer                               :: ileftw,iworkw,ilow,mflagw
  double precision                      :: basisu(ku),basisv(kv),basisw(kw)
  double precision                      :: worku(2*ku),workv(2*kv),workw(2*kw)

  val(:,:,:) = 0.0
  ilou = 1
  ilov = 1
  ilow = 1
  
  do ii=1,n
     do jj=1,m
        ! U
        call INTRV(tu,nctlu+ku,u(ii,jj),ilou,ileftu,mflagu)
        if (mflagu == 1) then
           ileftu = ileftu-ku
        end if
        call BSPVN(tu,ku,ku,1,u(ii,jj),ileftu,basisu,worku,iworku)
        istartu = ileftu-ku
     
        ! V
        call INTRV(tv,nctlv+kv,v(ii,jj),ilov,ileftv,mflagv)
        if (mflagv == 1) then
           ileftv = ileftv-kv
        end if
        call BSPVN(tv,kv,kv,1,v(ii,jj),ileftv,basisv,workv,iworkv)
        istartv = ileftv-kv
        
        ! W
        call INTRV(tw,nctlw+kw,w(ii,jj),ilow,ileftw,mflagw)
        if (mflagw == 1) then
           ileftw = ileftw-kw
        end if
        call BSPVN(tw,kw,kw,1,w(ii,jj),ileftw,basisw,workw,iworkw)
        istartw = ileftw-kw
        
        do i=1,ku
           do j=1,kv
              do k=1,kw
                 do idim=1,ndim
                    val(ii,jj,idim) = val(ii,jj,idim) + basisu(i)*basisv(j)*basisw(k)*coef(istartu+i,istartv+j,istartw+k,idim)
                 end do
              end do
           end do
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
  
  implicit none
  ! Input
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,ndim,n,m,l
  double precision, intent(in)          :: u(n,m,l),v(n,m,l),w(n,m,l)
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  double precision, intent(in)          :: coef(nctlu,nctlv,nctlw,ndim)

  ! Output
  double precision, intent(out)         :: val(n,m,l,ndim)
! Working
  integer                               :: idim,istartu,istartv,istartw,i,j,k,ii,jj,kk
  integer                               :: ileftu,iworku,ilou,mflagu
  integer                               :: ileftv,iworkv,ilov,mflagv
  integer                               :: ileftw,iworkw,ilow,mflagw
  double precision                      :: basisu(ku),basisv(kv),basisw(kw)
  double precision                      :: worku(2*ku),workv(2*kv),workw(2*kw)

  val(:,:,:,:) = 0.0
  ilou = 1
  ilov = 1
  ilow = 1

  do ii=1,n
     do jj=1,m
        do kk=1,l
           ! U
           call INTRV(tu,nctlu+ku,u(ii,jj,kk),ilou,ileftu,mflagu)
           if (mflagu == 1) then
              ileftu = ileftu-ku
           end if
           call BSPVN(tu,ku,ku,1,u(ii,jj,kk),ileftu,basisu,worku,iworku)
           istartu = ileftu-ku
           
           ! V
           call INTRV(tv,nctlv+kv,v(ii,jj,kk),ilov,ileftv,mflagv)
           if (mflagv == 1) then
              ileftv = ileftv-kv
           end if
           call BSPVN(tv,kv,kv,1,v(ii,jj,kk),ileftv,basisv,workv,iworkv)
           istartv = ileftv-kv
           
           ! W
           call INTRV(tw,nctlw+kw,w(ii,jj,kk),ilow,ileftw,mflagw)
           if (mflagw == 1) then
              ileftw = ileftw-kw
           end if
           call BSPVN(tw,kw,kw,1,w(ii,jj,kk),ileftw,basisw,workw,iworkw)
           istartw = ileftw-kw
           
           do i=1,ku
              do j=1,kv
                 do k=1,kw
                    do idim=1,ndim
                       val(ii,jj,kk,idim) = val(ii,jj,kk,idim) + basisu(i)*basisv(j)*basisw(k)*coef(istartu+i,istartv+j,istartw+k,idim)
                    end do
                 end do
              end do
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
