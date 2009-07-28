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
  double precision, intent(in)     :: s(nref)
  double precision, intent(in)     :: t(nref+2)
  double precision, intent(in)     :: x(nref,3)
  double precision, intent(in)     :: rot(nref,3)
  double precision, intent(in)     :: scale(nref)

  ! Link Data
  double precision, intent(in)     :: s_pos(nx,ny)
  double precision, intent(in)     :: links(nx,ny,3)

  ! Output
  double precision, intent(out)    :: coef(nx,ny,3)

  ! Local
  double precision                 :: local_s
  double precision                 :: matx(3,3)
  double precision                 :: maty(3,3)
  double precision                 :: matz(3,3)
  double precision                 :: rot_mat(3,3)
  double precision                 :: inv_rot_mat(3,3)

  double precision                 :: X_base(3)
  double precision                 :: angle,current_s,current_scale,det

  double precision bvalu
  double precision                :: work(3*2)

  print *,'in getcoef'
  i=1
  inbv = 1
  if (dir .eq. 1) then  ! Ref axis along V
     do j = 1,Ny
        print *,'along v'
        print *,'i,j,s_pos:',i,j,s_pos(1,j)
        current_s = s_pos(1,j)
        print *,'current_s:',current_s
        print *,'t:',t
        print *,'scale:',scale
        print *,'nref:',nref
        
        current_scale = bvalu(t,scale,nref,2,0,current_s,inbv,work)        
        ! Now we need the rotation Matrix
        print *,'done scale'
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
        print *,'done inverse'
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

        current_s = s_pos(i,1)
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
