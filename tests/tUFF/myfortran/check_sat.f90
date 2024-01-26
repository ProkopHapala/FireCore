subroutine check_sat ( bond_orders, changed )
  use conf, only: natoms, nbonds, neighbors, connectivity, valence_target
  implicit none
  integer, intent(inout) :: bond_orders(nbonds)
  logical, intent(out) :: changed
  integer :: ia, ib, j, nset, val
  logical :: found
  
  changed = .false.
  
  do ia = 1, natoms
     ! skip already set atoms
     found = .false.
     do j = 1, neighbors(ia)
        call find_bond ( ia, connectivity(ia,j), ib )
        if ( bond_orders(ib) < 0 ) found = .true.
     end do
     if ( .not. found ) cycle
     ! look for atoms that misses only one bond (or only single bonds)
     nset = 0
     val = 0
     do j = 1, neighbors(ia)
        call find_bond ( ia, connectivity(ia,j), ib )
        if ( bond_orders(ib) > 0 ) then
           nset = nset + 1
           val = val + bond_orders(ib)
        end if
     end do
     ! if found, set it
     if ( nset == neighbors(ia)-1 .or. neighbors(ia)-nset == valence_target(ia)-val ) then
        changed = .true.
        if ( nset == neighbors(ia)-1 ) then
           do j = 1, neighbors(ia)
              call find_bond ( ia, connectivity(ia,j), ib )
              if ( bond_orders(ib) < 0 ) then
                 bond_orders(ib) = valence_target(ia) - val
                 exit
              end if
           end do
        else if ( neighbors(ia)-nset == valence_target(ia)-val ) then
          do j = 1, neighbors(ia)
              call find_bond ( ia, connectivity(ia,j), ib )
              if ( bond_orders(ib) < 0 ) bond_orders(ib) = 1
           end do
        end if
     end if
  end do

  return
end subroutine check_sat
