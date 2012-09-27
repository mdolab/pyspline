subroutine getBasisPtSurface(u,v,tu,tv,ku,kv,vals,col_ind,istart,l_index,nctlu,nctlv,nnz)
  use precision
  implicit none

  ! Input
  integer         , intent(in)          :: ku,kv,nctlu,nctlv,nnz,istart
  real(kind=realType), intent(in)          :: u,v
  real(kind=realType), intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
  integer         , intent(in)          :: l_index(Nctlv,Nctlu)
  real(kind=realType), intent(inout)       :: vals(nnz)
  integer         , intent(inout)       :: col_ind(nnz)
  ! Working
  real(kind=realType)                      :: basisu(ku),basisv(kv)
  integer                               :: ileftu,mflagu,ilou
  integer                               :: ileftv,mflagv,ilov
  integer                               :: i,j,ii,jj,counter,start

  ilou = 1
  ilov = 1
  ! Get u interval
  call intrv(tu,nctlu+ku,u,ilou,ileftu,mflagu)
  if (mflagu == 0) then
     call basis(tu,nctlu,ku,u,ileftu,basisu)
  else if (mflagu == 1) then
     ileftu = nctlu
     basisu(:) = 0.0
     basisu(ku) = 1.0
  end if

  ! Get v interval
  call intrv(tv,nctlv+kv,v,ilov,ileftv,mflagv)
  if (mflagv == 0) then
     call basis(tv,nctlv,kv,v,ileftv,basisv)
  else if (mflagv == 1) then
     ileftv = nctlv
     basisv(:) = 0.0
     basisv(kv) = 1.0
  end if
  counter = 0
  do ii=1,ku
     do jj = 1,kv
        ! Get the local row/col for this surface
        start = istart + counter + 1
        col_ind(start) = l_index(ileftv-kv+jj,ileftu-ku+ii)
        vals(start) = basisu(ii)*basisv(jj)
        counter = counter + 1
    end do
 end do
end subroutine getBasisPtSurface

subroutine getBasisPtVolume(u,v,w,tu,tv,tw,ku,kv,kw,vals,col_ind,istart,l_index,nctlu,nctlv,nctlw,nnz)
  use precision
  implicit none
  ! Input
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,nnz,istart
  real(kind=realType), intent(in)          :: u,v,w
  real(kind=realType), intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  integer         , intent(in)          :: l_index(Nctlw,Nctlv,Nctlu)
  real(kind=realType), intent(inout)       :: vals(nnz)
  integer         , intent(inout)       :: col_ind(nnz)
  ! Working
  real(kind=realType)                      :: basisu(ku),basisv(kv),basisw(kw)
  integer                               :: ileftu,mflagu,ilou
  integer                               :: ileftv,mflagv,ilov
  integer                               :: ileftw,mflagw,ilow

  integer                               :: i,j,ii,jj,kk,counter,start

  ilou = 1
  ilov = 1
  ilow = 1
  ! Get u interval

  call intrv(tu,nctlu+ku,u,ilou,ileftu,mflagu)
  if (mflagu == 0) then
     call basis(tu,nctlu,ku,u,ileftu,basisu)
  else if (mflagu == 1) then
     ileftu = nctlu
     basisu(:) = 0.0
     basisu(ku) = 1.0
  end if

  ! Get v interval
  call intrv(tv,nctlv+kv,v,ilov,ileftv,mflagv)
  if (mflagv == 0) then
     call basis(tv,nctlv,kv,v,ileftv,basisv)
  else if (mflagv == 1) then
     ileftv = nctlv
     basisv(:) = 0.0
     basisv(kv) = 1.0
  end if

  ! Get w interval
  call intrv(tw,nctlw+kw,w,ilow,ileftw,mflagw)
  if (mflagw == 0) then
     call basis(tw,nctlw,kw,w,ileftw,basisw)
  else if (mflagw == 1) then
     ileftw = nctlw
     basisw(:) = 0.0
     basisw(kw) = 1.0
  end if

  counter = 0
  do ii=1,ku
     do jj = 1,kv
        do kk = 1,kw
           ! Get the local row/col for this surface
           start = istart + counter + 1
           
           col_ind(start) = l_index(ileftw-kw+kk,ileftv-kv+jj,ileftu-ku+ii)

           vals(start) = basisu(ii)*basisv(jj)*basisw(kk)
           counter = counter + 1
        end do
     end do
  end do
end subroutine getBasisPtVolume

