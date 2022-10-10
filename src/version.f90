
subroutine version(ver)
  ! Automatically generated version function
  !
  ! Desricption of Arguments
  ! Output
  ! ver    - string, contains git tag and commit hash
  use precision
  implicit none
  character(len=256), intent(out) :: ver
  
  ver = GITVERSION

end subroutine version
