subroutine assign_C1
  use conf, only: natoms, set_atom, neighbors, ufftypes, tipo!, connectivity, set_bond, bond_orders, bond_orders_int
  implicit none
  integer :: ia!, ib, j, ja
  
  do ia = 1, natoms
     if ( set_atom(ia) ) cycle
     if ( tipo(ia) == 'C' .and. neighbors(ia) == 2 ) then
        ufftypes(ia) = 'C_1'
        set_atom(ia) = .true.
!        ! set both bonds as cumulene, maybe it helps to finalize saturation, anyway it will be fixed by special_case_cumulene
!        do j = 1, 2
!           ja = connectivity(ia,j)
!           call find_bond ( ia, ja, ib )
!           bond_orders(ib) = 2.d0
!           bond_orders_int(ib) = 2
!           set_bond(ib) = .true.
!        end do
     end if
  end do
     
  return
end subroutine assign_C1
