subroutine assign_bonds
  use conf, only: nbonds, bonds, ufftypes, tol, nbondtypes, bondtypesmap, bond_orders
  use pars, only: bond_pars
  use uff, only: Zstar_I
  implicit none
  integer :: i, j, ib, m
  character(5) :: ati, atj
  logical :: found
  character(5) :: tmp(nbonds,2)
  real(8) :: tmp2(nbonds)

  ! calculate number of different bond types and store type arrays
  allocate ( bondtypesmap(nbonds), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for BONDTYPESMAP'
  nbondtypes = 0
  do ib = 1, nbonds
     found = .false.
     do m = 1, nbondtypes
        if ( ( ( ufftypes(bonds(ib,1)) == tmp(m,1) .and. ufftypes(bonds(ib,2)) == tmp(m,2) ) .or. &
             ( ufftypes(bonds(ib,1)) == tmp(m,2) .and. ufftypes(bonds(ib,2)) == tmp(m,1) ) ) .and. &
             abs(bond_orders(ib)-tmp2(m)) < tol ) then
           found = .true.
           bondtypesmap(ib) = m
           exit
        end if
     end do
     if ( .not. found ) then
        nbondtypes = nbondtypes + 1
        tmp(nbondtypes,1:2) = ufftypes(bonds(ib,1:2))
        tmp2(nbondtypes) = bond_orders(ib)
        bondtypesmap(ib) = nbondtypes
     end if
  end do
  
  ! compute pars
  allocate ( bond_pars(nbondtypes,2), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for BOND_PARS'
  do m = 1, nbondtypes
     ati = tmp(m,1)
     atj = tmp(m,2)
     call find_I ( ati, i )
     call find_I ( atj, j )
     call calc_r_IJ ( i, j, tmp2(m), bond_pars(m,2) )
     bond_pars(m,1) = 0.5d0 * 664.12d0 * Zstar_I(i) * Zstar_I(j) / bond_pars(m,2)**3
     write(*,'(a,i4,2(1x,a))') 'bond type ', m, trim(ati), trim(atj)
  end do
  
  return
end subroutine assign_bonds
