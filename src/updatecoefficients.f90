subroutine updatecoefficients(type,s_pos,links,coef,indicies,s,t,x,&
     rotx,roty,rotz,scale,rot_type,nref,ns,ncoef)
!*** DESCRIPTION
!
!     Written by Gaetan Kenway
! 
!     The function takes the important data from the reference axis,
!     as well as the links for the surface and calculates the new
!     positions of the control points. This is a complex-only function
!     as it is only used for complex-step calculations
!
!     Description of Arguments:
!     Input:
!     type    - integer: 0 if stored links are used, 1 if 'full' type is used 
!               and links must be recomputed before rotation can be applied
!     s_pos: Real, vector, size Ns,  The s position connection along direction of ref axis
!     links: Real, matrix, size Ns by 3:  The local distance connections 
!     coef    - Real, List of displaced coefficients, size ncoef by 3
!     Ref Axis Data:
!          s  - Real, vector of length nref - the basis for the ref axis
!          t  - Real, vector of length nref + 2: Knots for ref axis
!          x  - Real, matrix of length nref by 3: Control points for the 
!               spatial part of ref axis
!          rot - Real, matrix of length nref by 3: Control points for the 
!               rotation part of ref axis
!          scale - Real, vector or length nref: Control points for the 
!                  scale data on ref axis
!
!     Output:
!     coef    - Real, List of displaced coefficients, size ncoef by 3
!
  implicit none

  real, PARAMETER :: pi = 3.14159265358979

  integer ns,i,idim,nref,inbv,ncoef
  integer         , intent(in)     :: type,rot_type
  integer         , intent(in)     :: indicies(ns)
  ! Ref axis:
  double precision, intent(in)     :: s(nref)
  double precision, intent(in)     :: t(nref+2)
  double precision, intent(in)     :: x(nref,3)
  double precision, intent(in)     :: rotx(nref),roty(nref),rotz(nref)
  double precision, intent(in)     :: scale(nref)

  ! Link Data
  double precision, intent(in)     :: s_pos(ns)
  double precision, intent(inout)  :: links(ns,3)

  ! Output
  double precision, intent(inout)  :: coef(ncoef,3)

  ! Local
  double precision                 :: local_s
  double precision                 :: matx(3,3)
  double precision                 :: maty(3,3)
  double precision                 :: matz(3,3)
  double precision                 :: rot_mat(3,3)
  double precision                 :: inv_rot_mat(3,3)
  double precision                 :: X_base(3)
  double precision                 :: angle,current_s,current_scale,det

  do i = 1,ns
     current_s = s_pos(i)
     call eval_curve(current_s,t,2,scale,nref,1,current_scale)

     ! X Rot:
     call eval_curve(current_s,t,2,rotx,nref,1,angle)
     call xrot(angle*pi/180,matx)

     ! Y Rot:
     call eval_curve(current_s,t,2,roty,nref,1,angle)
       call yrot(angle*pi/180,maty)

     ! Z Rot:
     call eval_curve(current_s,t,2,rotz,nref,1,angle)
     call zrot(angle*pi/180,matz)
     !Select the rotation order
     if (rot_type == 1) then
        rot_mat = transpose(matmul(matz,matmul(maty,matx)))
     else if (rot_type == 2) then
        rot_mat = transpose(matmul(maty,matmul(matz,matx)))
     else if (rot_type == 3) then
        rot_mat = transpose(matmul(matx,matmul(matz,maty)))          
     else if (rot_type == 4) then
        rot_mat = transpose(matmul(matz,matmul(matx,maty)))          
     else if (rot_type == 5) then
        rot_mat = transpose(matmul(maty,matmul(matx,matz)))          
     else if (rot_type == 6) then
        rot_mat = transpose(matmul(matx,matmul(maty,matz)))          
     end if
   
     call eval_curve(current_s,t,2,x,nref,3,X_base)
     if (type == 1) then ! We have a full connection and need to recompute links
        links(i,:) = coef(indicies(i)+1,:) - X_base
     end if
     
     coef(indicies(i)+1,:) = X_base + matmul(rot_mat,links(i,:))*current_scale
  end do
end subroutine updatecoefficients


subroutine xrot(theta,rotmat)
  !form rotation maxtrix for x-axis
  implicit none
  double precision ::  theta
  double precision ::  rotmat(3,3)

  rotmat(1,1) = 1
  rotmat(1,2) = 0
  rotmat(1,3) = 0
  
  rotmat(2,1) = 0
  rotmat(2,2) = cos(theta)
  rotmat(2,3) = -sin(theta)
  
  rotmat(3,1) = 0
  rotmat(3,2) = sin(theta)
  rotmat(3,3) = cos(theta)
  
end subroutine xrot

subroutine yrot(theta,rotmat)
  !form rotation maxtrix for y-axis
  implicit none
  double precision ::  theta
  double precision ::  rotmat(3,3)
  
  rotmat(1,1) = cos(theta)
  rotmat(1,2) = 0
  rotmat(1,3) = -sin(theta)
  
  rotmat(2,1) = 0
  rotmat(2,2) = 1
  rotmat(2,3) = 0
  
  rotmat(3,1) = sin(theta)
  rotmat(3,2) = 0
  rotmat(3,3) = cos(theta)
  
end subroutine yrot

subroutine zrot(theta,rotmat)
  !form rotation maxtrix for z-axis
  implicit none
  double precision ::  theta
  double precision ::  rotmat(3,3)
  
  rotmat(1,1) = cos(theta)
  rotmat(1,2) = -sin(theta)
  rotmat(1,3) = 0
  
  rotmat(2,1) = sin(theta)
  rotmat(2,2) = cos(theta)
  rotmat(2,3) = 0
  
  rotmat(3,1) = 0
  rotmat(3,2) = 0
  rotmat(3,3) = 1
  
end subroutine zrot
