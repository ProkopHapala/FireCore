subroutine conjugation
  use conf, only: natoms, neighbors, connectivity, tipo, ufftypes, bond_orders
  implicit none
  integer :: ia, j, ja, ib
  
  do ia = 1, natoms
     if ( ( tipo(ia) == 'N' .and. neighbors(ia) == 3 ) .or. &
          ( tipo(ia) == 'O' .and. neighbors(ia) == 2 ) ) then
        do j = 1, neighbors(ia)
           ja = connectivity(ia,j)
           !if ( ufftypes(ja)(3:3) == 'R' .or. ufftypes(ja) == 'C_2' .or. ufftypes(ja) == 'N_2' ) then
           if ( ufftypes(ja)(3:3) == 'R' ) then
              if ( tipo(ia) == 'N' ) ufftypes(ia) = 'N_R'
              if ( tipo(ia) == 'O' ) ufftypes(ia) = 'O_R'
              call find_bond ( ia, ja, ib )
              bond_orders(ib) = 1.5d0
           end if
        end do
     end if
  end do

  return
end subroutine conjugation
