subroutine quality_volume(X,nx,ny,nz,quality)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract quality_volume evaluates the quality of the volumes inside X
  !              (SCLAR VERSION)
  !     Description of Arguments
  !     Input
  !     X       - Real, Coordinates: nx,ny,nz,3
  !     nx      - Real, Number of coordinates in x
  !     ny      - Real, Number of coordinates in y
  !     nz      - Real, Number of coordinates in z
  !
  !     Ouput 
  !     quality - list of quality of volumes, size (nx-1)*(ny-1)*(nz-1)
  
  implicit none
  ! Input
  integer         , intent(in)          :: nx,ny,nz
  double precision, intent(in)          :: X(nx,ny,nz,3)
  ! Output
  double precision, intent(out)         :: quality((nx-1)*(ny-1)*(nz-1))
  
  ! Working
  integer :: counter ,indices(8,3),i,j,k,ii
  double precision :: points(8,3)
  counter = 0

  do i=1,nx-1
     do j=1,ny-1
        do k=1,nz-1
           counter = counter + 1
           
           call hexa_index(i,j,k,indices)
           
           do ii=1,8
              points(ii,:) = X(indices(ii,1),indices(ii,2),indices(ii,3),:)
           end do
           call skew_hexa(points,quality(counter))
        end do
     end do
  end do
end subroutine quality_volume

subroutine quality_volume_deriv(X,nx,ny,nz,offset,localIndex,vals,col_ind,nQuality)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract quality_volume evaluates the quality of the volumes inside X
  !              (SCLAR VERSION)
  !     Description of Arguments
  !     Input
  !     X       - Real, Coordinates: nx,ny,nz,3
  !     nx      - Real, Number of coordinates in x
  !     ny      - Real, Number of coordinates in y
  !     nz      - Real, Number of coordinates in z
  !     offset  - Integer: The row to start filling up array
  !     localIndex - Integer the local to global mapping for X, size nx,ny,nz
  !     vals    - Real, CSR array of values; size: 24*nQuality
  !     col_ind - integer, CSR array of col indices; size: 24*nQuality
  !     nQuality- Integer, the total number of coef volumes in pyBlock object
  !     Ouput 
  !     vals    - As above
  !     col_ind - As above
  
  implicit none
  ! Input
  integer         , intent(in)          :: nx,ny,nz
  double precision, intent(in)          :: X(nx,ny,nz,3)
  integer         , intent(in)          :: offset,nQuality
  integer         , intent(in)          :: localIndex(nx,ny,nz)
  double precision, intent(inout)       :: vals(24*nQuality)
  integer         , intent(inout)       :: col_ind(24*nQuality)

  ! Working
  integer :: counter ,indices(8,3),i,j,k,ii,jj
  double precision :: points(8,3),pointsb(8,3)
  double precision :: quality,qualityb

  counter = 0

   do i=1,nx-1
      do j=1,ny-1
         do k=1,nz-1

           call hexa_index(i,j,k,indices)
           do ii=1,8
              points(ii,:) = X(indices(ii,1),indices(ii,2),indices(ii,3),:)
           end do
           qualityb = 1.0
           call skew_hexa_b(points,pointsb,quality,qualityb)
        
           ! pointsb has the derivative we're after
           do ii=1,8
              do jj=1,3
                 counter = counter + 1
                 vals(offset+counter) = pointsb(ii,jj)
                 col_ind(offset+counter) = 3*localIndex(indices(ii,1),indices(ii,2),indices(ii,3)) + jj -1
              end do
           end do

        end do
     end do
  end do
end subroutine quality_volume_deriv

subroutine verify_quality_volume_deriv(X,nx,ny,nz)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract quality_volume evaluates the quality of the volumes inside X
  !              (SCLAR VERSION)
  !     Description of Arguments
  !     Input
  !     X       - Real, Coordinates: nx,ny,nz,3
  !     nx      - Real, Number of coordinates in x
  !     ny      - Real, Number of coordinates in y
  !     nz      - Real, Number of coordinates in z
  
  implicit none
  ! Input
  integer         , intent(in)          :: nx,ny,nz
  double precision, intent(in)          :: X(nx,ny,nz,3)
  
  ! Working
  integer :: counter ,indices(8,3),i,j,k,ii,jj
  double precision :: points(8,3),pointsb(8,3)
  double precision :: quality,qualityb,h,qualityd(8,3),quality0

  counter = 0
  h = 1.0e-5
   do i=1,nx-1
      do j=1,ny-1
         do k=1,nz-1
            counter = counter + 1
          
           call hexa_index(i,j,k,indices)
           do ii=1,8
              points(ii,:) = X(indices(ii,1),indices(ii,2),indices(ii,3),:)
           end do
           qualityb = 1.0
           call skew_hexa_b(points,pointsb,quality,qualityb)
        
           ! pointsb has the derivative we're after
           call skew_hexa(points,quality0)
           do ii=1,8
              do jj=1,3
                 points(ii,jj) = points(ii,jj) + h
                 call skew_hexa(points,quality)
                 qualityd(ii,jj) = (quality-quality0)/h
                 points(ii,jj) = points(ii,jj) - h
              end do
           end do

           print *,'Derivative for i,j,k:',i,j,k
           do ii=1,8
              do jj=1,3
                 print *,ii,jj,qualityd(ii,jj),pointsb(ii,jj)
              end do
           end  do

        end do
     end do
  end do
end subroutine verify_quality_volume_deriv
  
subroutine skew_hexa(points,skew)

  implicit none
  double precision  :: points(8,3),skew
  double precision  :: scalar_triple_product
  skew = 1.0
  
  skew = skew*scalar_triple_product(points(2,:)-points(1,:),points(3,:)-points(1,:),points(5,:)-points(1,:))
  skew = skew*scalar_triple_product(points(4,:)-points(2,:),points(1,:)-points(2,:),points(6,:)-points(2,:))
  skew = skew*scalar_triple_product(points(1,:)-points(3,:),points(4,:)-points(3,:),points(7,:)-points(3,:))
  skew = skew*scalar_triple_product(points(3,:)-points(4,:),points(2,:)-points(4,:),points(8,:)-points(4,:))
  skew = skew*scalar_triple_product(points(7,:)-points(5,:),points(6,:)-points(5,:),points(1,:)-points(5,:))
  skew = skew*scalar_triple_product(points(5,:)-points(6,:),points(8,:)-points(6,:),points(2,:)-points(6,:))
  skew = skew*scalar_triple_product(points(8,:)-points(7,:),points(5,:)-points(7,:),points(3,:)-points(7,:))
  skew = skew*scalar_triple_product(points(6,:)-points(8,:),points(7,:)-points(8,:),points(4,:)-points(8,:))
end subroutine skew_hexa

function scalar_triple_product(a,b,c)
  implicit none 
  double precision :: a(3),b(3),c(3),scalar_triple_product
  ! Compute the scalar triple product of a*(b x c)
  ! This is the determinate of:

  !                   [a1 a2 a3] [a b c]
  !(b x c) * a  = det [b1 b2 b3] [d e f]
  !                   [c1 c2 c3] [g h i]
  !                         aei          +        bfg     +        cdh     −     afh       −     bdi       − ceg.
  a = a/(sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3)))
  b = b/(sqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3)))
  c = c/(sqrt(c(1)*c(1)+c(2)*c(2)+c(3)*c(3)))
  scalar_triple_product = a(1)*b(2)*c(3) + a(2)*b(3)*c(1) + a(3)*b(1)*c(2) - a(1)*b(3)*c(2) - a(2)*b(1)*c(3) - a(3)*b(2)*c(1)

end function scalar_triple_product
  
subroutine hexa_index(i,j,k,indices)

  integer :: i,j,k,indices(8,3)

  indices(1,:) = (/i,j,k/)
  indices(2,:) = (/i+1,j,k/)
  indices(3,:) = (/i,j+1,k/)
  indices(4,:) = (/i+1,j+1,k/)
  indices(5,:) = (/i,j,k+1/)
  indices(6,:) = (/i+1,j,k+1/)
  indices(7,:) = (/i,j+1,k+1/)
  indices(8,:) = (/i+1,j+1,k+1/)
end subroutine hexa_index
