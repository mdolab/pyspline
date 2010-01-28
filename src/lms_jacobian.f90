module lms_jacobian

  implicit none
  public   :: setup_jacobian,kill_jacobian,Aprod1,Aprod2
  
  ! DYNAMIC WORKSPACE DEFINED HERE.
  ! They are allocated in lsqrtest and used by Aprod1, Aprod2.
  
 
  integer         , allocatable ::  col_ind(:),row_ptr(:)
  double precision, allocatable ::  vals(:)

contains

  subroutine setup_jacobian(nrow,ncol,nnz_row)
    integer          ::  nrow,ncol,nnz,nnz_row

    ! Setup the lms_jacobain object
    nnz = nrow*nnz_row
    allocate(col_ind(nnz))
    allocate(row_ptr(nrow))
    allocate(vals(nnz))
    
  end subroutine setup_jacobian

  subroutine kill_jacobian()
    deallocate(col_ind,row_ptr,vals)
  end subroutine kill_jacobian

  subroutine Aprod1(m,n,x,y)

    integer,          intent(in)    :: m,n
    double precision, intent(in)    :: x(n)
    double precision, intent(inout) :: y(m)

    !-------------------------------------------------------------------
    ! Aprod1 computes y = y + A*x without altering x,
    ! where A is a test matrix of the form  A = Y*D*Z,
    ! and the matrices D, Y, Z are represented by
    ! the allocatable vectors d, hy, hz in this module.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------

    ! Do a sparse Matrix-Vector Produc
    
    
  end subroutine Aprod1

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Aprod2(m,n,x,y)

    integer,          intent(in)    :: m,n
    double precision, intent(in)    :: x(n)
    double precision, intent(inout) :: y(m)

    !-------------------------------------------------------------------
    ! Aprod2 computes x = x + A'*y without altering y,
    ! where A is a test matrix of the form  A = Y*D*Z,
    ! and the matrices D, Y, Z are represented by
    ! the allocatable vectors d, hy, hz in this module.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------
    
    ! Do a sparse Matrix-Vector product

  end subroutine Aprod2

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module lms_jacobian
