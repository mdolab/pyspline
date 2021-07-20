
subroutine curve_jacobian_wrap(s, sd, t, k, nctl, n, nd, vals, row_ptr, col_ind)
  use precision
  implicit none
  ! Input
  integer              , intent(in)      :: k, nctl, n, nd
  real(kind=realType)  , intent(in)      :: t(nctl+k), s(n), sd(nd)
  ! Output
  real(kind=realType)  , intent(inout)   :: vals((n+nd)*k)
  integer              , intent(inout)   :: row_ptr(n+nd+1)
  integer              , intent(inout)   :: col_ind((n+nd)*k)
  real(kind=realType)                   :: basisu(k), basisud(k,k)
  integer                               :: i, j, counter, ileft

  counter = 1
  do i=1, n ! Do the values first 
     call findSpan(s(i), k, t, nctl, ileft)
     call basis(t, nctl, k, s(i), ileft, basisu)

     row_ptr(i) = counter-1
     do j=1, k
        col_ind(counter) = ileft-k+j-1
        vals(counter) = basisu(j)
        counter = counter + 1
     end do
  end do
  do i=1, nd ! Do the derivatives next
     call findSpan(sd(i), k, t, nctl, ileft)
     call derivBasis(t, nctl, k, sd(i), ileft, 1, basisud)

     row_ptr(i+n) = counter-1
     do j=1, k
        col_ind(counter) = ileft-k+j-1
        vals(counter) = basisud(2, j)
        counter = counter + 1
     end do
  end do
  row_ptr(n+nd+1) = counter-1
end subroutine curve_jacobian_wrap

subroutine constr_jac(A_val, A_row_ptr, A_col_ind, B_val, B_row_ptr, B_col_ind, C_val, C_row_ptr, C_col_ind, &
     Am, An, Cm, Annz, Bnnz, Cnnz, J_val, J_col_ind, J_row_ptr)
  use precision
  implicit none
  ! Input
  integer         , intent(in)   :: Am, An, Cm, Annz, Bnnz, Cnnz
  real(kind=realType), intent(in)   :: A_val(Annz), B_val(Bnnz), C_val(Cnnz)
  integer         , intent(in)   :: A_col_ind(Annz), B_col_ind(Bnnz), C_col_ind(Cnnz)
  integer         , intent(in)   :: A_row_ptr(Am+1), B_row_ptr(Am+1), C_row_ptr(Cm+1)

  ! Output
  real(kind=realType), intent(out)  :: J_val(Annz+Bnnz+Cnnz)
  integer         , intent(out)  :: J_col_ind(Annz+Bnnz+Cnnz)
  integer         , intent(out)  :: J_row_ptr(Am+Cm+1)

  ! Local 
  integer                        :: i, j, counter
  ! This functions assembes the following CSR matrix:
  ! J = [A    B]
  !     [C    0]

  ! Now assmeble the full jacobain
  counter = 1
  J_row_ptr(1) = 0
  do i =1, Am
     do j=1, A_row_ptr(i+1)-A_row_ptr(i)
        J_val(counter) =     A_val(A_row_ptr(i)+j)
        J_col_ind(counter) = A_col_ind(A_row_ptr(i)+j)
        counter = counter + 1
     end do
     do j=1, B_row_ptr(i+1)-B_row_ptr(i)
        J_val(counter) = B_val(B_row_ptr(i)+j)
        J_col_ind(counter) = B_col_ind(B_row_ptr(i)+j) + An
        counter = counter + 1
     end do
     J_row_ptr(i+1) = counter - 1
  end do
  do i =1, Cm
     do j=1, C_row_ptr(i+1)-C_row_ptr(i)
        J_val(counter) = C_val(C_row_ptr(i)+j)
        J_col_ind(counter) = C_col_ind(C_row_ptr(i)+j)
        counter = counter + 1
     end do
     J_row_ptr(i+1+am) = counter - 1
  end do
end subroutine constr_jac

subroutine poly_length(X, n, ndim, length)
  ! Compute the length of the spatial polygon
  use precision
  implicit none

  !Input
  integer             , intent(in)    :: n, ndim
  real(kind=realType) , intent(in)    :: X(ndim, n)
  
  ! Ouput
  real(kind=realType), intent(out)   :: length

  ! Working
  integer                         :: i, idim
  real(kind=realType)             :: dist

  length = 0.0
  do i=1, n-1
     dist = 0.0
     do idim=1, ndim
        dist = dist + (X(idim, i)-X(idim, i+1))**2
     end do
     length = length + sqrt(dist)
  end do
  
end subroutine poly_length

subroutine curve_para_corr(t, k, s, coef, nctl, ndim, length, n, X)

  ! Do Hoschek parameter correction
  use precision
  implicit none
  ! Input/Output
  integer              , intent(in)      :: k, nctl, ndim, n
  real(kind=realType)  , intent(in)      :: t(nctl+k)
  real(kind=realType)  , intent(inout)   :: s(n)
  real(kind=realType)  , intent(in)      :: coef(ndim, nctl)
  real(kind=realType)  , intent(in)      :: X(ndim, n)
  real(kind=realType)  , intent(in)      :: length
  ! Working
  integer                               :: i, j, max_inner_iter
  real(kind=realType)                   :: D(ndim), D2(ndim)
  real(kind=realType)                   :: val(ndim), deriv(ndim)
  real(kind=realType)                   :: c, s_tilde

  max_inner_iter = 10
  do i=2, n-1
     call eval_curve(s(i), t, k, coef, nctl, ndim, 1, val)
     call eval_curve_deriv(s(i), t, k, coef, nctl, ndim, deriv)
          
     D = X(:, i)-val
     c = dot_product(D, deriv)/NORM2(deriv)

     inner_loop: do j=1, max_inner_iter

        s_tilde = s(i)+ c*(t(nctl+k)-t(1))/length
        call eval_curve(s_tilde, t, k, coef, nctl, ndim, 1, val)
        D2 = X(:, i)-val
        if (NORM2(D) .ge. NORM2(D2)) then
           s(i) = s_tilde
           exit inner_loop
        else
           c = c*0.5
        end if
     end do inner_loop
  end do

end subroutine curve_para_corr

function compute_rms_curve(t, k, s, coef, nctl, ndim, n, X)
  ! Compute the rms
  use precision
  implicit none
  ! Input/Output
  integer           , intent(in)      :: k, nctl, ndim, n
  real(kind=realType)  , intent(in)      :: t(k+nctl)
  real(kind=realType)  , intent(in)      :: s(n)
  real(kind=realType)  , intent(in)      :: coef(ndim, nctl)
  real(kind=realType)  , intent(in)      :: X(ndim, n)
  real(kind=realType)                   :: compute_rms_curve 

  ! Working
  integer                            :: i, idim
  real(kind=realType)                :: val(ndim), D(ndim)
  
  compute_rms_curve = 0.0
  do i=1, n
     call eval_curve(s(i), t, k, coef, nctl, ndim, 1, val)
     D = val-X(:, i)
     do idim=1, ndim
        compute_rms_curve = compute_rms_curve + D(idim)**2
     end do
  end do
  compute_rms_curve = sqrt(compute_rms_curve/n)
end function compute_rms_curve
