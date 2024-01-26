subroutine check_aroma ( natoms, atoms, aromatic )
  use conf, only: tol, bond_orders, bond_orders_int, neighbors, tipo
  implicit none
  integer, intent(in) :: natoms, atoms(0:natoms+1)
  logical, intent(out) :: aromatic
  integer :: ia, ib1, ib2, npi, npi_save

  aromatic = .false.
  
  ! check that there is no sp3 carbons in the ring
  do ia = 1, natoms 
     if ( tipo(atoms(ia)) == 'C' .and. neighbors(atoms(ia)) /= 3 ) then
        aromatic = .false.
        return
     end if
  end do
  
  npi = 0

  ! count heteroatoms in the ring that can donate an electron pair
  do ia = 1, natoms
     if ( ( tipo(atoms(ia)) == 'N' .and. neighbors(atoms(ia)) == 3 ) .or. &
          ( tipo(atoms(ia)) == 'O' .and. neighbors(atoms(ia)) == 2 ) ) npi = npi + 2
  end do

  npi_save = npi
  
  ! count pi electrons in the ring considering double bonds in the ring only
  do ia = 1, natoms
     call find_bond ( atoms(ia), atoms(ia+1), ib1 )
     call find_bond ( atoms(ia-1), atoms(ia), ib2 )
     if ( bond_orders_int(ib1) == 2 .or. bond_orders_int(ib2) == 2 ) npi = npi + 1
  end do

  if ( npi == 6 ) then
     aromatic = .true.
     return
  end if

  npi = npi_save

  ! count pi electrons in the ring considering also fused rings
  do ia = 1, natoms
     call find_bond ( atoms(ia), atoms(ia+1), ib1 )
     call find_bond ( atoms(ia-1), atoms(ia), ib2 )
     if ( bond_orders_int(ib1) == 2 .or. abs(bond_orders(ib1)-1.5d0) < tol .or. &
          bond_orders_int(ib2) == 2 .or. abs(bond_orders(ib2)-1.5d0) < tol ) npi = npi + 1
  end do

  if ( npi == 6 ) aromatic = .true.
  
  return
end subroutine check_aroma
