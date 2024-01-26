subroutine special_case_nitro
  use conf, only: natoms, bond_orders, set_atom, set_bond, &
       neighbors, connectivity, ufftypes, bond_orders, tipo
  implicit none
  integer :: ia, ib, j, ja, n
  
  do ia = 1, natoms
     !if ( set_atom(ia) ) cycle
     if ( tipo(ia) == 'N' .and. neighbors(ia) == 3 ) then
        n = 0
        do j = 1, 3
           ja = connectivity(ia,j)
           if ( tipo(ja) == 'O' ) n = n + 1
        end do
        if ( n == 2 ) then
           ufftypes(ia) = 'N_R'
           set_atom(ia) = .true.
           do j = 1, 3
              ja = connectivity(ia,j)
              if ( tipo(ja) == 'O' ) then
                 ufftypes(ja) = 'O_R'
                 set_atom(ja) = .true.
                 call find_bond ( ia, ja, ib )
                 bond_orders(ib) = 1.5d0
                 set_bond(ib) = .true.
              else
                 call find_bond ( ia, ja, ib )
                 bond_orders(ib) = 1.d0
                 set_bond(ib) = .true.
              end if
           end do
        end if
     end if
  end do

  return
end subroutine special_case_nitro
