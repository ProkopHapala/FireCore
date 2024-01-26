subroutine assign_C2
  use conf, only: natoms, set_atom, neighbors, ufftypes, tipo
  implicit none
  integer :: ia
  
  do ia = 1, natoms
     if ( set_atom(ia) ) cycle
     if ( tipo(ia) == 'C' .and. neighbors(ia) == 3 ) then
        ufftypes(ia) = 'C_2'
        set_atom(ia) = .true.
     end if
  end do
     
  return
end subroutine assign_C2
