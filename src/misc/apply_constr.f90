subroutine apply_constr(t,n,k,Kstiff,constr,ndim,nconstr,Knew,f)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: apply_constr applies the explict constraints to the
  !     stiffness matrix and generates the rhs vector
  !
  !     Description of Arguments
  !     Input
  !     t       - Real, Vector of knots, size n+k
  !     n       - Integer, number of control points
  !     k       - Integer, spline order
  !     Kstiff  - Real, Matrix of size nxn, Computed Stiffness
  !     const   - Real, Matrix of size nconstr,ndim+2 - Constraint matrix. (i,1) is the s position
  !               constr(i,2) is the x-displacement and constr(i,3) is the y-displacement
  !     ndim    -Integer, the sptatial dimension, usually 2 or 3
  !     nconstr -Integer, size of the const matrix
  !     Output
  !     Knew    - Real, Matrix of size nxn, New stiffness matrix
  !     f       - Real, Vector of size n, Forcing rhs vector

  !    Input
  double precision,   intent(in) :: t(n+k)
  integer ,           intent(in) :: n,k,nconstr,ndim
  double precision,   intent(in) :: Kstiff(n,n)
  double precision,   intent(in) :: constr(nconstr,ndim+2)

  !    Output
  double precision,  intent(out) :: Knew(n,n)
  double precision,  intent(out) :: f(n,ndim)
  !    Working
  integer i,j,iconstr
  double precision               :: B(n)
  external get_basis
  print *, 'Hello'
  ! Copy the stiffness passed in to the new matrix
  do i=1,n
     do j =i,n
        Knew(i,j) = Kstiff(i,j)
     end do
  end do

  do iconstr=1,nconstr
     ! Get the basis vectors for this constraint
     do i=1,n
        B(i) = get_basis(t,n,k,i,constr(iconstr,1),0)
     end do

     do i=1,n
        Knew(i,:) = Knew(i,:) + B(i)*B*constr(iconstr,2+ndim)
     end do

     do idim =1,ndim
        f(:,idim) = f(:,idim) + constr(iconstr,2+ndim)*constr(iconstr,1+idim)*B
     end do
  end do
end subroutine apply_constr
