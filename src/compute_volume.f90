subroutine volume_jacobian_wrap(u,v,w,tu,tv,tw,ku,kv,kw,nctlu,nctlv,nctlw,nu,nv,nw,vals,row_ptr,col_ind)

  implicit none
  ! Input
  integer         , intent(in)          :: ku,kv,kw,nctlu,nctlv,nctlw,nu,nv,nw
  double precision, intent(in)          :: u(nu,nv,nw),v(nu,nv,nw)
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv),tw(nctlw+kw)
  double precision, intent(out)         :: vals(nu*nv*nw*ku*kv*kw)
  integer         , intent(out)         :: col_ind(nu*nv*nw*ku*kv*kw),row_ptr(nu*nv*nw+1)
  ! Working
  double precision                      :: basisu(ku),basisv(kv),basisw(kw)
  integer                               :: ilou,ileftu,mflagu,worku(2*ku),iworku
  integer                               :: ilov,ileftv,mflagv,workv(2*kv),iworkv
  integer                               :: ilow,ileftw,mflagw,workw(2*kw),iworkw
  integer                               :: i,j,k,ii,jj,kk,counter,c1,c2

  ilou = 1
  ilov = 1
  ilow = 1
  counter = 1

  do i=1,nu
     do j = 1,nv
        do k=1,kw
           ! U
           call INTRV(tu,nctlu+ku,u(i,j,k),ilou,ileftu,mflagu)
           if (mflagu == 1) then
              ileftu = ileftu-ku
           end if
           call BSPVN(tu,ku,ku,1,u(i,j,k),ileftu,basisu,worku,iworku)

           ! V
           call INTRV(tv,nctlv+kv,v(i,j,k),ilov,ileftv,mflagv)
           if (mflagv == 1) then
              ileftv = ileftv-kv
           end if
           call BSPVN(tv,kv,kv,1,v(i,j,k),ileftv,basisv,workv,iworkv)

           ! W
           call INTRV(tw,nctlw+kw,w(i,j,k),ilow,ileftw,mflagw)
           if (mflagw == 1) then
              ileftw = ileftw-kw
           end if
           call BSPVN(tw,kw,kw,1,w(i,j,k),ileftw,basisw,workw,iworkw)
        
           row_ptr( (i-1)*nv + j  ) = counter-1
           do ii=1,ku
              c1 = (ileftu-ku+ii-1)*Nctlv*Nctlw
              do jj = 1,kv
                 c2 = (ileftv-kv+jj-1)*Nctlw
                 do kk=1,kw
                    col_ind(counter) = c1 + c2 + (ileftw-kw+kk-1)
                    vals(counter) = basisu(ii)*basisv(jj)*basisw(kk)
                    counter = counter + 1
                 end do
              end do
           end do
        end do
     end do
  end do
  row_ptr(nu*nv*nw+1) = counter-1
  
end subroutine volume_jacobian_wrap
