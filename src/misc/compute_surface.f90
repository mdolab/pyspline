subroutine surface_jacobian_linear(u,v,tu,tv,ku,kv,nctlu,nctlv,nu,nv,rows,cols,vals)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract surface_jacobian_linear computes the non-zero values in 
  !     in the jacobain matrix which are returned to be set in a petsc4py 
  !     sparse matrix
  !
  !     Description of Arguments
  !     Input
  !     u       - Real, u coordinates, size: nu x nv
  !     v       - Real, v coordinates, size: nu x nv
  !     tu      - Real,Knot vector in u. Length nctlu+ku
  !     tv      - Real,Knot vector in v. Length nctlv+kv
  !     ku      - Integer, order of B-spline in u
  !     kv      - Integer, order of B-spline in v
  !     coef    - Real,Array of B-spline coefficients  Size (nctlu,nctlv,ndim)
  !     nctlu   - Integer,Number of control points in u
  !     nctlv   - Integer,Number of control points in v
  !     nu      - Integer, Number of data points in u
  !     nv      - Integer, Number of data points in v
  !     ndim    - Integer, Spatial Dimension
  !
  !     Ouput 
  !     rows    - Integer, Row index, length n
  !     cols    - Integer, Col index, length n
  !     vals    - Real   , Values, length n

  implicit none
  ! Input
  integer         , intent(in)          :: ku,kv,nctlu,nctlv,nu,nv
  double precision, intent(in)          :: u(nu,nv),v(nu,nv)
  double precision, intent(in)          :: tu(nctlu+ku),tv(nctlv+kv)
  ! Output
  integer          ,intent(out)         :: rows(nu*nv*ku*kv),cols(nu*nv*ku*kv)
  double precision, intent(out)         :: vals(nu*nv*ku*kv)

  ! Working
  double precision                      :: vniku(ku),worku(2*ku)
  integer                               :: ilou,ileftu,mflagu

  double precision                      :: vnikv(kv),workv(2*kv)
  integer                               :: ilov,ileftv,mflagv

  integer                               :: i,j,ii,jj,counter,iwork

  ilou = 1
  ilov = 1
  counter = 1
  do i=1,nu
     do j = 1,nv
        ! Get u interval
        call intrv(tu,nctlu+ku,u(i,j),ilou,ileftu,mflagu)
        if (mflagu == 0) then
           call bspvn(tu,ku,ku,1,u(i,j),ileftu,vniku,worku,iwork)
        else if (mflagu == 1) then
           ileftu = nctlu
           vniku(:) = 0.0
           vniku(ku) = 1.0
        end if

        ! Get v interval
        call intrv(tv,nctlv+kv,v(i,j),ilov,ileftv,mflagv)
        if (mflagv == 0) then
           call bspvn(tv,kv,kv,1,v(i,j),ileftv,vnikv,workv,iwork)
        else if (mflagv == 1) then
           ileftv = nctlv
           vnikv(:) = 0.0
           vnikv(kv) = 1.0
        end if

        do ii=1,ku
           do jj = 1,kv
              rows(counter) = (i-1)*nv + j -1
              cols(counter) = (ileftu-ku+ii-1)*Nctlu + (ileftv-kv+jj-1)
              vals(counter) = vniku(ii)*vnikv(jj)
              counter = counter + 1
           end do
        end do
     end do
  end do
end subroutine surface_jacobian_linear

subroutine surface_para_corr(tu,tv,ku,kv,u,v,coef,nctlu,nctlv,ndim,nu,nv,X,rms)

  ! Do Hoschek parameter correction
  implicit none
  ! Input/Output
  double precision  ,intent(in)      :: tu(ku+nctlu),tv(kv+nctlv)
  double precision  ,intent(inout)   :: u(nu,nv),v(nu,nv)
  double precision  ,intent(in)      :: coef(nctlu,nctlv,ndim)
  integer           ,intent(in)      :: ku,kv,nctlu,nctlv,ndim,nu,nv
  double precision  ,intent(in)      :: X(nu,nv,ndim)
  double precision  ,intent(out)     :: rms
  ! Working
  integer                            :: i,j,jj,max_inner_iter
  double precision                   :: lengthu,lengthv
  double precision                   :: D(ndim),D2(ndim),Dnorm,D2norm
  double precision                   :: val(ndim),deriv(2,ndim)
  double precision                   :: delta_c,delta_d,u_tilde,v_tilde
  integer                            :: adj_u,adj_v
  !Functions
  double precision                   :: norm,poly_length

  max_inner_iter = 10
  rms = 0.0
  do i=1,nu
     do j = 1,nv

        lengthu = poly_length(X(:,j,:),nu,ndim)
        lengthv = poly_length(X(i,:,:),nv,ndim)

        adj_u = 1 ! Adjust them by default
        adj_v = 1

        if (i==1 .or. i==nu) then
           adj_v = 0
        end if

        if (j==1 .or. j == nv) then
           adj_u = 0
        end if

        call eval_surface(u(i,j),v(i,j),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
        call eval_surface_deriv(u(i,j),v(i,j),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,deriv)

        D = X(i,j,:)-val

        if (adj_u == 1) then
           delta_c = dot_product(D,deriv(1,:))/norm(deriv(1,:),ndim)

           inner_loop1: do jj=1,max_inner_iter
              u_tilde = u(i,j)+ delta_c*(tu(nctlu+ku)-tu(1))/lengthu
              call eval_surface(u_tilde,v(i,j),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
              D2 = X(i,j,:)-val
              if (norm(D,ndim) .ge. norm(D2,ndim)) then
                 u(i,j) = u_tilde
                 exit inner_loop1
              else
                 delta_c = delta_c*0.5
              end if
           end do inner_loop1
        end if

        if (adj_v == 1) then
           delta_d = dot_product(D,deriv(2,:))/norm(deriv(2,:),ndim)

           inner_loop2: do jj=1,max_inner_iter
              v_tilde = v(i,j)+ delta_d*(tv(nctlv+kv)-tv(1))/lengthv
              call eval_surface(u(i,j),v_tilde,tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
              D2 = X(i,j,:)-val
              if (norm(D,ndim) .ge. norm(D2,ndim)) then
                 v(i,j) = v_tilde
                 exit inner_loop2
              else
                 delta_d = delta_d*0.5
              end if
           end do inner_loop2
        end if

     end do
  end do

  ! Lets redo the full RMS
  rms = 0.0
  do i=1,nu
     do j=1,nv
        call eval_surface(u(i,j),v(i,j),tu,tv,ku,kv,coef,nctlu,nctlv,ndim,val)
        D = X(i,j,:)-val
        rms = rms + dot_product(D,D)
     end do
  end do
  rms = sqrt(rms/(nu*nv))
  

end subroutine surface_para_corr

