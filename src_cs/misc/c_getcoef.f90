subroutine getcoef(dir,s,t,x,rot,scale,s_pos,links,coef,nref,nx,ny)

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
!     dir     - Integer: 0 or 1: Determines along which direction the
!               reference axis is attached. 0 for along u, 1 along v
!     Ref Axis Data:
!          s  - Real, vector of length nref - the basis for the ref axis
!          t  - Real, vector of length nref + 2: Knots for ref axis
!          x  - Real, matrix of length nref by 3: Control points for the 
!               spatial part of ref axis
!          rot - Real, matrix of length nref by 3: Control points for the 
!               rotation part of ref axis
!          scale - Real, vector or length nref: Control points for the 
!                  scale data on ref axis
!     s_pos: Real, vector, size Ns,  The s position connection along direction of ref axis
!     links: Real, matrix, size Ns by 3:  The local distance connections 
!     


!     Output:
!     coef    - Real, Matrix of displaced coefficients, size nx by ny
!
  implicit none

  real, PARAMETER :: pi = 3.14159265358979

  integer nx,ny,i,j,idim,nref,inbv
  integer   , intent(in)     :: dir
  ! Ref axis:
  complex*16, intent(in)     :: s(nref)
  complex*16, intent(in)     :: t(nref+2)
  complex*16, intent(in)     :: x(nref,3)
  complex*16, intent(in)     :: rot(nref,3)
  complex*16, intent(in)     :: scale(nref)

  ! Link Data
  complex*16, intent(in)     :: s_pos(nx,ny)
  complex*16, intent(in)     :: links(nx,ny,3)

  ! Output
  complex*16, intent(out)    :: coef(nx,ny,3)

  ! Local
  complex*16                 :: local_s
  complex*16                 :: matx(3,3)
  complex*16                 :: maty(3,3)
  complex*16                 :: matz(3,3)
  complex*16                 :: rot_mat(3,3)
  complex*16                 :: inv_rot_mat(3,3)

  complex*16                 :: X_base(3)
  complex*16                 :: angle,current_s,current_scale,det

  complex*16 bvalu
  complex*16                :: work(3*2)

  print *,'in getcoef'
  
  if (dir .eq. 1) then  ! Ref axis along V
     do j = 1,Ny
        !print *,'i,j,s_pos:',i,j,s_pos(1,j)
        current_s = s_pos(1,j)
        current_scale = bvalu(t,scale,nref,2,0,current_s,inbv,work)        
        ! Now we need the rotation Matrix

        ! Python Command:
        !inv(dot(self._roty(self.rotys(s)), dot(self._rotx(self.rotxs(s)),self._rotz(self.rotzs(s)))))
        
        ! X Rot:
        angle = bvalu(t,rot(:,1),nref,2,0,current_s,inbv,work)
        call xrot(angle*pi/180,matx)

        ! Y Rot:
        angle = bvalu(t,rot(:,2),nref,2,0,current_s,inbv,work)
        call yrot(angle*pi/180,maty)

        ! Z Rot:
        angle = bvalu(t,rot(:,3),nref,2,0,current_s,inbv,work)
        call zrot(angle*pi/180,matz)

        rot_mat = matmul(maty,matmul(matx,matz))


        rot_mat = matmul(maty,matmul(matx,matz))
        
        ! DO the inverse of the 3x3 explicitly
        det = rot_mat(3,3)*( rot_mat(1,1)*rot_mat(2,2) - rot_mat(2,1)*rot_mat(1,2) ) &
        - rot_mat(3,2)*( rot_mat(1,1)*rot_mat(2,3) - rot_mat(2,1)*rot_mat(1,3) ) &
        + rot_mat(3,1)*( rot_mat(1,2)*rot_mat(2,3) - rot_mat(1,3)*rot_mat(2,2) )
    
        inv_rot_mat(1,1) =   (rot_mat(2,2)*rot_mat(3,3) - rot_mat(2,3)*rot_mat(3,2))/det
        inv_rot_mat(1,2) = - (rot_mat(1,2)*rot_mat(3,3) - rot_mat(1,3)*rot_mat(3,2))/det
        inv_rot_mat(1,3) =   (rot_mat(1,2)*rot_mat(2,3) - rot_mat(1,3)*rot_mat(2,2))/det
    
        inv_rot_mat(2,1) = - (rot_mat(2,1)*rot_mat(3,3) - rot_mat(2,3)*rot_mat(3,1))/det
        inv_rot_mat(2,2) =   (rot_mat(1,1)*rot_mat(3,3) - rot_mat(1,3)*rot_mat(3,1))/det
        inv_rot_mat(2,3) = - (rot_mat(1,1)*rot_mat(2,3) - rot_mat(1,3)*rot_mat(2,1))/det
        
        inv_rot_mat(3,1) =   (rot_mat(2,1)*rot_mat(3,2) - rot_mat(2,2)*rot_mat(3,1))/det
        inv_rot_mat(3,2) = - (rot_mat(1,1)*rot_mat(3,2) - rot_mat(1,2)*rot_mat(3,1))/det
        inv_rot_mat(3,3) =   (rot_mat(1,1)*rot_mat(2,2) - rot_mat(1,2)*rot_mat(2,1))/det
        
        do idim =1,3
           X_base(idim) = bvalu(t,x(:,idim),nref,2,0,current_s,inbv,work)
        end do
        
        do i =1,Nx
           coef(i,j,:) = X_base + matmul(inv_rot_mat,links(i,j,:))*current_scale
        end do
     end do 

  else
     i=1
     j=1

     do i = 1,Nx

        current_s = s_pos(i,j)
        current_scale = bvalu(t,scale,nref,2,0,current_s,inbv,work)        
        ! Now we need the rotation Matrix

        ! Python Command:
        !inv(dot(self._roty(self.rotys(s)), dot(self._rotx(self.rotxs(s)),self._rotz(self.rotzs(s)))))
        
        ! X Rot:
        angle = bvalu(t,rot(:,1),nref,2,0,current_s,inbv,work)
        call xrot(angle*pi/180,matx)

        ! Y Rot:
        angle = bvalu(t,rot(:,2),nref,2,0,current_s,inbv,work)
        call yrot(angle*pi/180,maty)

        ! Z Rot:
        angle = bvalu(t,rot(:,3),nref,2,0,current_s,inbv,work)
        call zrot(angle*pi/180,matz)

        rot_mat = matmul(maty,matmul(matx,matz))
        
        ! DO the inverse of the 3x3 explicitly
        det = rot_mat(3,3)*( rot_mat(1,1)*rot_mat(2,2) - rot_mat(2,1)*rot_mat(1,2) ) &
        - rot_mat(3,2)*( rot_mat(1,1)*rot_mat(2,3) - rot_mat(2,1)*rot_mat(1,3) ) &
        + rot_mat(3,1)*( rot_mat(1,2)*rot_mat(2,3) - rot_mat(1,3)*rot_mat(2,2) )
    
        inv_rot_mat(1,1) =   (rot_mat(2,2)*rot_mat(3,3) - rot_mat(2,3)*rot_mat(3,2))/det
        inv_rot_mat(1,2) = - (rot_mat(1,2)*rot_mat(3,3) - rot_mat(1,3)*rot_mat(3,2))/det
        inv_rot_mat(1,3) =   (rot_mat(1,2)*rot_mat(2,3) - rot_mat(1,3)*rot_mat(2,2))/det
    
        inv_rot_mat(2,1) = - (rot_mat(2,1)*rot_mat(3,3) - rot_mat(2,3)*rot_mat(3,1))/det
        inv_rot_mat(2,2) =   (rot_mat(1,1)*rot_mat(3,3) - rot_mat(1,3)*rot_mat(3,1))/det
        inv_rot_mat(2,3) = - (rot_mat(1,1)*rot_mat(2,3) - rot_mat(1,3)*rot_mat(2,1))/det
        
        inv_rot_mat(3,1) =   (rot_mat(2,1)*rot_mat(3,2) - rot_mat(2,2)*rot_mat(3,1))/det
        inv_rot_mat(3,2) = - (rot_mat(1,1)*rot_mat(3,2) - rot_mat(1,2)*rot_mat(3,1))/det
        inv_rot_mat(3,3) =   (rot_mat(1,1)*rot_mat(2,2) - rot_mat(1,2)*rot_mat(2,1))/det
        
        do idim =1,3
           X_base(idim) = bvalu(t,x(:,1),nref,2,0,current_s,inbv,work)
        end do
        
        
        do j =1,Ny
           coef(i,j,:) = X_base + matmul(inv_rot_mat,links(i,j,:))*current_scale
        end do
     end do 


  end if


end subroutine getcoef


subroutine xrot(theta,rotmat)
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
  
end subroutine xrot

subroutine yrot(theta,rotmat)
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
  
end subroutine yrot

subroutine zrot(theta,rotmat)
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
  
end subroutine zrot
