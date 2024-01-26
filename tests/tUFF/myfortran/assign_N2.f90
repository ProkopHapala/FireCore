subroutine assign_N2
  use conf, only: natoms, set_atom, neighbors, ufftypes, tipo
  implicit none
  integer :: ia
  
  do ia = 1, natoms
     if ( set_atom(ia) ) cycle
     if ( tipo(ia) == 'N' .and. neighbors(ia) == 2 ) then
        ufftypes(ia) = 'N_2'
        set_atom(ia) = .true.
     end if
  end do
     
  return
end subroutine assign_N2
