subroutine simple_rule
  use conf, only: tol, natoms, ufftypes, neighbors, connectivity, bond_orders, &
       set_atom, set_bond, tipo
  implicit none
  integer :: ia, ib, ja, j, n
  
  
  ! simplest-est rule: atom with two sp2 neighbors (or heteroatoms) is resonant
  do ia = 1, natoms
     if ( ( tipo(ia) == 'C' .and. neighbors(ia) == 3 ) .or. & ! C_2
          ( tipo(ia) == 'N' .and. neighbors(ia) > 1 ) .or. &  ! N_3 or N_2
          ( tipo(ia) == 'O' .and. neighbors(ia) == 2 ) ) then ! O_3
        n = 0
        do j = 1, neighbors(ia)
           ja = connectivity(ia,j)
           if ( ( tipo(ja) == 'C' .and. neighbors(ja) == 3 ) .or. & ! C_2
                ( tipo(ja) == 'N' .and. neighbors(ja) > 1 ) .or. &  ! N_3 or N_2
                ( tipo(ja) == 'O' .and. neighbors(ja) == 2 ) ) n = n + 1
        end do
        if ( n > 1 ) then
           ! set atoms
           if ( set_atom(ia) ) then
              if ( ufftypes(ia)(3:3) /= 'R' ) write(*,*) 'WARNING SIMPLERULE: atom ', tipo(ia), ia, &
                   ' would be set to resonant but it has already a type of ', ufftypes(ia)
              cycle
           end if
           if ( tipo(ia) == 'C' ) then
              ufftypes(ia) = 'C_R'
              set_atom(ia) = .true.
           end if
           if ( tipo(ia) == 'N' ) then
              ufftypes(ia) = 'N_R'
              set_atom(ia) = .true.
           end if
           if ( tipo(ia) == 'O' ) then
              ufftypes(ia) = 'O_R'
              set_atom(ia) = .true.
           end if
           ! set bonds
           do j = 1, neighbors(ia)
              ja = connectivity(ia,j)
              if ( ( tipo(ja) == 'C' .and. neighbors(ja) == 3 ) .or. & ! C_2
                   ( tipo(ja) == 'N' .and. neighbors(ja) > 1 ) .or. &  ! N_3 or N_2
                   ( tipo(ja) == 'O' .and. neighbors(ja) == 2 ) ) then
                 call find_bond ( ia, ja, ib )
                 if ( set_bond(ib) ) then
                    if ( abs(bond_orders(ib)-1.5d0) > tol ) write(*,*) 'WARNING SIMPLERULE: bond ', ib, ' between atoms ', &
                         tipo(ia), ia, tipo(ja), ja, ' would be set to 1.5 but it has already a bond order of ', bond_orders(ib)
                    cycle
                 end if
                 bond_orders(ib) = 1.5d0
                 set_bond(ib) = .true.
              end if
           end do
        end if
     end if
  end do

  return
end subroutine simple_rule
