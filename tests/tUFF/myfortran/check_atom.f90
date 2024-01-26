subroutine check_atom ( ia, bond_orders, exclude )
  use conf, only: nbonds, neighbors, connectivity
  implicit none
  integer, intent(in) :: ia, bond_orders(nbonds)
  logical, intent(out) :: exclude
  integer :: ib, j
  
  exclude = .true.

  do j = 1, neighbors(ia)
     call find_bond ( ia, connectivity(ia,j), ib )
     if ( bond_orders(ib) < 0 ) then
        exclude = .false.
        return
     end if
  end do

  return
end subroutine check_atom
