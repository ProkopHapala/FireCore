subroutine special_case_amide
  use conf, only: natoms, bond_orders, set_atom, set_bond, neighbors, connectivity, &
       ufftypes, bond_orders, bond_orders_int, tipo
  implicit none
  integer :: ia, ib, j, ja
  logical :: found

  do ia = 1, natoms
     !if ( set_atom(ia) ) cycle
     if ( tipo(ia) == 'C' .and. neighbors(ia) == 3 ) then
        ! look for a carbonyl oxygen
        found = .false.
        do j = 1, neighbors(ia)
           ja = connectivity(ia,j)
           if ( tipo(ja) == 'O' .and. neighbors(ja) == 1 ) found = .true.
        end do
        if ( .not. found ) cycle
        ! look for an amino nitrogen
        found = .false.
        do j = 1, neighbors(ia)
           ja = connectivity(ia,j)
           if ( tipo(ja) == 'N' .and. neighbors(ja) == 3 ) found = .true.
        end do
        if ( .not. found ) cycle
        ufftypes(ia) = 'C_R'
        set_atom(ia) = .true.
        do j = 1, neighbors(ia)
           ja = connectivity(ia,j)
           if ( tipo(ja) == 'O' .and. neighbors(ja) == 1 ) then
              ufftypes(ja) = 'O_R'
              set_atom(ja) = .true.
              call find_bond ( ia, ja, ib )
              bond_orders(ib) = 1.5d0
              bond_orders_int(ib) = 2
              set_bond(ib) = .true.
           else if ( tipo(ja) == 'N' .and. neighbors(ja) == 3 ) then
              ufftypes(ja) = 'N_R'
              set_atom(ja) = .true.
              call find_bond ( ia, ja, ib )
              !bond_orders(ib) = 1.5d0 ! manually change resonant bonds between C_R and N_R from 1.5 to 1.41
              bond_orders(ib) = 1.41d0 ! manually change resonant bonds between C_R and N_R from 1.5 to 1.41
              bond_orders_int(ib) = 1
              set_bond(ib) = .true.
              ! commented out because of amide group in aromatic rings
!           else
!              call find_bond ( ia, ja, ib )
!              bond_orders(ib) = 1.d0
!              bond_orders_int(ib) = 1
!              set_bond(ib) = .true.
           end if
        end do
     end if
  end do

  return
end subroutine special_case_amide
