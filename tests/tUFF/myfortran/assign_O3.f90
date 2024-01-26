subroutine assign_O3
  use conf, only: natoms, bond_orders, set_atom, set_bond, neighbors, connectivity, &
       ufftypes, bond_orders, bond_orders_int, tipo
  implicit none
  integer :: ia, ib, j, ja
  
  do ia = 1, natoms
     if ( set_atom(ia) ) cycle
     if ( tipo(ia) == 'O' .and. neighbors(ia) == 2 ) then
        ufftypes(ia) = 'O_3'
        set_atom(ia) = .true.
        do j = 1, 2
           ja = connectivity(ia,j)
           call find_bond ( ia, ja, ib )
           bond_orders(ib) = 1.d0
           bond_orders_int(ib) = 1
           set_bond(ib) = .true.
        end do
     end if
  end do
  
  return
end subroutine assign_O3
