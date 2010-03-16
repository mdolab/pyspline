SUBROUTINE basis(t,nctl,k,x,ileft,B)

  implicit none
  integer nctl,k,ileft,l,j,ipj,jp1ml,jp1
  double precision :: t(nctl+k),B(k),vm,vmprev,x,diff1(k),diff2(k)
    
  B(:) = 0.0
  B(1) = 1.0

  do j=1,k-1
     ipj = ileft + j
     diff1(j) = t(ipj)-x
     diff2(j) = x-t(ileft-j+1)
     vmprev = 0.0
     jp1 = j+1
     do l=1,j
        jp1ml = jp1-l
        vm = B(l)/(diff1(l)+diff2(jp1ml))
        B(l) = vm*diff1(l) + vmprev
        vmprev = vm*diff2(jp1ml)
     end do
     B(jp1) = vmprev
  end do
end SUBROUTINE basis
