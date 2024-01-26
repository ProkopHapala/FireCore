subroutine assign_O2
  use conf, only: natoms, bond_orders, set_atom, set_bond, neighbors, connectivity, &
       ufftypes, bond_orders, bond_orders_int, tipo
  implicit none
  integer :: ia, ib, ja
  
  do ia = 1, natoms
     if ( set_atom(ia) ) cycle
     if ( tipo(ia) == 'O' .and. neighbors(ia) == 1 ) then
        ufftypes(ia) = 'O_2'
        set_atom(ia) = .true.
        ja = connectivity(ia,1)
        call find_bond ( ia, ja, ib )
        bond_orders(ib) = 2.d0
        bond_orders_int(ib) = 2
        set_bond(ib) = .true.
        if ( tipo(ja) == 'C' .and. .not. set_atom(ja) .and. neighbors(ja) == 3 ) then
           ufftypes(ja) = 'C_2'
           set_atom(ja) = .true.
        end if
     end if
  end do
     
  return
end subroutine assign_O2
