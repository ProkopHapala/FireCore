subroutine find_bond ( ia1, ia2, ib )
  use conf, only: nbonds, bonds
  implicit none
  integer, intent(in) :: ia1, ia2
  integer, intent(out) :: ib
  integer :: i
  
  ib = -1

  do i = 1, nbonds
     if ( ( bonds(i,1) == ia1 .and. bonds(i,2) == ia2 ) .or. ( bonds(i,1) == ia2 .and. bonds(i,2) == ia1 ) ) then
        ib = i
        return
     end if
  end do

  write(0,*) 'ERROR FINDBOND: bond not found ', ia1, ia2
  write(0,*) 'STOP'
  stop
  
  return
end subroutine find_bond
