subroutine trivial_cases
  use conf, only: natoms, bond_orders, set_atom, set_bond, neighbors, connectivity, &
       ufftypes, bond_orders, bond_orders_int, tipo
  implicit none
  integer :: ia, ib, j, ja, ka
  
  ! NB: all bond_orders are set to -1.0, set_atom and set_bond to false
  
  ! set fixed atoms and bonds (i.e. atoms that are already saturated with sigma skeleton)
  do ia = 1, natoms
     if ( tipo(ia) == 'H' ) then ! all hydrogens
        if ( neighbors(ia) /= 1 ) then
           write(0,*) 'ERROR TRIVIALCASES: hydrogen ', ia, ' missing the bond'
           write(0,*) 'STOP'
           stop
        end if
        ufftypes(ia) = 'H_'
        set_atom(ia) = .true.
        call find_bond ( ia, connectivity(ia,1), ib )
        bond_orders(ib) = 1.d0
        bond_orders_int(ib) = 1
        set_bond(ib) = .true.
     else if ( tipo(ia) == 'C' .and. neighbors(ia) == 4 ) then ! carbons with 4 neighbors
        ufftypes(ia) = 'C_3'
        set_atom(ia) = .true.
        do j = 1, 4
           call find_bond ( ia, connectivity(ia,j), ib )
           bond_orders(ib) = 1.d0
           bond_orders_int(ib) = 1
           set_bond(ib) = .true.
        end do
     else if ( tipo(ia) == 'N' .and. neighbors(ia) == 3 ) then ! nitrogens with 3 neighbors
        ! not setting ufftype and final bond orders, as some may be resonant
        do j = 1, 3
           call find_bond ( ia, connectivity(ia,j), ib )
           bond_orders_int(ib) = 1
        end do
     else if ( tipo(ia) == 'N' .and. neighbors(ia) == 1 ) then ! nitrogens with 1 neighbor
        ufftypes(ia) = 'N_1'
        set_atom(ia) = .true.
        ja = connectivity(ia,1)
        call find_bond ( ia, ja, ib )
        bond_orders(ib) = 3.d0
        bond_orders_int(ib) = 3
        set_bond(ib) = .true.
        ! nitrile group
        if ( tipo(ja) == 'C' ) then
           ufftypes(ja) = 'C_1'
           set_atom(ja) = .true.
           do j = 1, 2
              ka = connectivity(ja,j)
              if ( ka /= ia ) then
                 call find_bond ( ja, ka, ib )
                 bond_orders(ib) = 1.d0
                 bond_orders_int(ib) = 1
                 set_bond(ib) = .true.
                 exit
              end if
           end do
        end if
     else if ( tipo(ia) == 'O' .and. neighbors(ia) == 2 ) then ! oxygens with 2 neighbors
        ! not setting ufftype and final bond orders, as they may be resonant
        do j = 1, 2
           call find_bond ( ia, connectivity(ia,j), ib )
           bond_orders_int(ib) = 1
        end do
     else if ( tipo(ia) == 'O' .and. neighbors(ia) == 1 ) then ! oxygens with 1 neighbor
        ! not setting ufftype and final bond order, as it may be resonant
        call find_bond ( ia, connectivity(ia,1), ib )
        bond_orders_int(ib) = 2
     end if
  end do

  return
end subroutine trivial_cases
