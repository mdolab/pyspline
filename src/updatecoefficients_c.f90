subroutine updatecoefficients_c(type,s_pos,links,coef,indicies,s,t,x,&
     rotx,roty,rotz,scale,rot_type,nref,ns,ncoef)
!*** DESCRIPTION
!
!     Written by Gaetan Kenway
! 
!     Complex verison for derivative calcs -- see updatecoefficients for description

  implicit none

  real, PARAMETER :: pi = 3.14159265358979

  integer ns,i,idim,nref,inbv,ncoef
  integer         , intent(in)     :: type,rot_type
  integer         , intent(in)     :: indicies(ns)
  ! Ref axis:
  double precision, intent(in)     :: s(nref)
  double precision, intent(in)     :: t(nref+2)
  complex*16, intent(in)     :: x(nref,3)
  complex*16, intent(in)     :: rotx(nref),roty(nref),rotz(nref)
  complex*16, intent(in)     :: scale(nref)

  ! Link Data
  double precision, intent(in)     :: s_pos(ns)
  double precision, intent(inout)  :: links(ns,3)

  ! Output
  complex*16, intent(inout)  :: coef(ncoef,3)

  ! Local
  double precision                 :: local_s
  complex*16                 :: matx(3,3)
  complex*16                 :: maty(3,3)
  complex*16                 :: matz(3,3)
  complex*16                 :: rot_mat(3,3)
  complex*16                 :: X_base(3)
  complex*16                 :: angle,current_s,current_scale

  do i = 1,ns
     current_s = s_pos(i)

     ! scale (the coefficients) are complex
     call eval_curve_c(current_s,t,2,scale,nref,1,current_scale)

     ! X Rot:
     call eval_curve_c(current_s,t,2,rotx,nref,1,angle)
     call xrot_c(angle*pi/180,matx)

     ! Y Rot:
     call eval_curve_c(current_s,t,2,roty,nref,1,angle)
       call yrot_c(angle*pi/180,maty)

     ! Z Rot:
     call eval_curve_c(current_s,t,2,rotz,nref,1,angle)
     call zrot_c(angle*pi/180,matz)
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
   
     call eval_curve_c(current_s,t,2,x,nref,3,X_base)
     if (type == 1) then ! We have a full connection and need to recompute links
        links(i,:) = coef(indicies(i)+1,:) - X_base
     end if
     
     coef(indicies(i)+1,:) = X_base + matmul(rot_mat,links(i,:))*current_scale
  end do
end subroutine updatecoefficients_c


subroutine xrot_c(theta,rotmat)
  !form rotation maxtrix for x-axis
  implicit none
  complex*16 ::  theta
  complex*16 ::  rotmat(3,3)

  rotmat(1,1) = 1
  rotmat(1,2) = 0
  rotmat(1,3) = 0
  
  rotmat(2,1) = 0
  rotmat(2,2) = cos(theta)
  rotmat(2,3) = -sin(theta)
  
  rotmat(3,1) = 0
  rotmat(3,2) = sin(theta)
  rotmat(3,3) = cos(theta)
  
end subroutine xrot_c

subroutine yrot_c(theta,rotmat)
  !form rotation maxtrix for y-axis
  implicit none
  complex*16 ::  theta
  complex*16 ::  rotmat(3,3)
  
  rotmat(1,1) = cos(theta)
  rotmat(1,2) = 0
  rotmat(1,3) = -sin(theta)
  
  rotmat(2,1) = 0
  rotmat(2,2) = 1
  rotmat(2,3) = 0
  
  rotmat(3,1) = sin(theta)
  rotmat(3,2) = 0
  rotmat(3,3) = cos(theta)
  
end subroutine yrot_c

subroutine zrot_c(theta,rotmat)
  !form rotation maxtrix for z-axis
  implicit none
  complex*16 ::  theta
  complex*16 ::  rotmat(3,3)
  
  rotmat(1,1) = cos(theta)
  rotmat(1,2) = -sin(theta)
  rotmat(1,3) = 0
  
  rotmat(2,1) = sin(theta)
  rotmat(2,2) = cos(theta)
  rotmat(2,3) = 0
  
  rotmat(3,1) = 0
  rotmat(3,2) = 0
  rotmat(3,3) = 1
  
end subroutine zrot_c
