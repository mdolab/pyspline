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
    allocate(row_ptr(nrow+1))
    allocate(vals(nnz))
    
  end subroutine setup_jacobian

  subroutine kill_jacobian()
    deallocate(col_ind,row_ptr,vals)
  end subroutine kill_jacobian

  subroutine Aprod1(m,n,x,y)
    use lsqrDataModule, only : dp
    integer,  intent(in)    :: m,n
    real(dp), intent(in)    :: x(n)
    real(dp), intent(inout) :: y(m)
    integer                         :: i,k1,k2
    !-------------------------------------------------------------------
    ! Aprod1 computes y = y + A*x without altering x,
    ! A is stored in compressed sparse row format
    !-------------------------------------------------------------------

    ! Do a sparse Matrix-Vector Produc
    do i=1,m
       k1 = row_ptr(i)
       k2 = row_ptr(i+1)-1
       Y(i) = Y(i) + dot_product(vals(k1:k2),x(col_ind(k1:k2)))
    end do
  end subroutine Aprod1

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Aprod2(m,n,x,y)
    use lsqrDataModule, only : dp
    integer,  intent(in)    :: m,n
    real(dp), intent(inout) :: x(n)
    real(dp), intent(in)    :: y(m)
    
    integer                         :: i,k
    !-------------------------------------------------------------------
    ! Aprod2 computes x = x + A'*y without altering y,
    ! A is stored in compress sparse row foramt
    !-------------------------------------------------------------------
    ! Do a sparse Matrix-Vector product
    
    do i=1,m
       do k=row_ptr(i),row_ptr(i+1)-1
          x(col_ind(k)) = x(col_ind(k)) + vals(k)*y(i)
       end do
    end do
    

  end subroutine Aprod2

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module lms_jacobian
