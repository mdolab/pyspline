subroutine getBasisPt(u,v,tu,tv,ku,kv,vals,col_ind,istart,l_index,nctlu,nctlv,nnz)

  implicit none
  ! Input
  integer         , intent(in)          :: ku,kv,nctlu,nctlv,nnz,istart
  double precision, intent(in)          :: u,v
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
  integer         , intent(in)          :: l_index(Nctlu,Nctlv)
  double precision, intent(inout)       :: vals(nnz)
  integer         , intent(inout)       :: col_ind(nnz)
  ! Working
  double precision                      :: vniku(ku),worku(2*ku)
  integer                               :: ilou,ileftu,mflagu,idim

  double precision                      :: vnikv(kv),workv(2*kv)
  integer                               :: ilov,ileftv,mflagv,start,g_index
  integer                               :: i,j,ii,jj,iwork,counter

  ilou = 1
  ilov = 1
  ! Get u interval
  call intrv(tu,nctlu+ku,u,ilou,ileftu,mflagu)
  if (mflagu == 0) then
     call bspvn(tu,ku,ku,1,u,ileftu,vniku,worku,iwork)
  else if (mflagu == 1) then
     ileftu = nctlu
     vniku(:) = 0.0
     vniku(ku) = 1.0
  end if

  ! Get v interval
  call intrv(tv,nctlv+kv,v,ilov,ileftv,mflagv)
  if (mflagv == 0) then
     call bspvn(tv,kv,kv,1,v,ileftv,vnikv,workv,iwork)
  else if (mflagv == 1) then
     ileftv = nctlv
     vnikv(:) = 0.0
     vnikv(kv) = 1.0
  end if
  counter = 0
  do idim = 0,2
     do ii=1,ku
        do jj = 1,kv
           ! Get the local row/col for this surface
           g_index = l_index(ileftu - ku + ii , ileftv - kv + jj )
           start = istart + counter + 1
           vals(start) = vniku(ii)*vnikv(jj)
           col_ind(start) = g_index*3 + idim
           counter = counter + 1
        end do
     end do
  end do
end subroutine getBasisPt

