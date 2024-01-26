subroutine find_I ( atI, I )
  use uff, only: nufftypes, sym
  implicit none
  character(5), intent(in) :: atI
  integer, intent(out) :: I
  integer :: j

  do j = 1, nufftypes
     if ( sym(j) == atI ) then
        I = j
        return
     end if
  end do
  
  write(0,*) 'ERROR FINDI: atom ', atI, ' not found'
  write(0,*) 'STOP'
  stop

end subroutine find_I
